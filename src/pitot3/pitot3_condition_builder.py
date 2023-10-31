#! /usr/bin/env python3

"""
Starting some condition builder messing...

Chris James (c.james4@uq.edu.au) - 02/01/23

"""

CONDITION_BUILDER_VERSION_STRING = '31-Oct-2023'

#----------------------------------------------------------------------------------------

import yaml
import numpy as np
import os
import sys
import itertools
from contextlib import redirect_stdout, redirect_stderr
from natsort import natsorted
import json
import zipfile
import shutil
import time

from pitot3 import run_pitot3
from pitot3_utils.pitot3_classes import pitot3_remote_run_creator, pitot3_pickle_output_file_loader, \
    pitot3_single_line_output_file_creator, pitot3_pickle_output_file_creator, pitot3_json_output_file_loader, Facility

#----------------------------------------------------------------------------------------

def run_pitot3_condition_builder(config_dict = {}, config_filename = None,
                                 pitot3_data_folder = '$PITOT3_DATA', default_config_yaml_filename = 'PITOT3_default_config.yaml',
                                 condition_builder_cleanup = False):

    print('-' * 60)
    print (f'Running PITOT3 Condition Builder Version: {CONDITION_BUILDER_VERSION_STRING}.')
    print('-' * 60)

    # --------------------------------------------------------------------------------
    # load the user's config

    if config_filename:

        user_config_yaml_filename = config_filename

        user_config_file = open(user_config_yaml_filename)
        user_config_data = yaml.load(user_config_file, Loader=yaml.FullLoader)

        for variable in user_config_data.keys():
            if user_config_data[variable] == 'None':
                user_config_data[variable] = None

    # --------------------------------------------------------------------------------
    # now merge the two config files from the input file and/or the config_dict
    # we put the file config and the dictionary input config both here, with the dictionary updated last,
    # so it will overwrite file values with dictionary ones if both user inputs are used.

    condition_builder_config_dict = {}
    if config_filename:
        condition_builder_config_dict.update(user_config_data)
    condition_builder_config_dict.update(config_dict)

    # --------------------------------------------------------------------------------
    # load the PITOT3 default config which loads everything which runs the program...

    default_config_file = open(os.path.expandvars(pitot3_data_folder + '/' + default_config_yaml_filename))
    default_config_data = yaml.load(default_config_file, Loader=yaml.FullLoader)

    # go through and any remove any None values which yaml has loaded as None
    for variable in default_config_data.keys():
        if default_config_data[variable] == 'None':
            default_config_data[variable] = None

    # we will remove anything from this which is used to set up the condition builder...
    base_config_dict = condition_builder_config_dict.copy()

    base_output_filename = condition_builder_config_dict['output_filename']

    if 'condition_builder_cleanup' in condition_builder_config_dict:
        condition_builder_cleanup = condition_builder_config_dict['condition_builder_cleanup']
    else:
        condition_builder_cleanup = False

    if condition_builder_cleanup:
        print('-' * 60)
        print("Removing any evidence of an old condition builder run with the same name from the folder.")
        print(f"(This will remove any file or folder related to the basename of {base_output_filename}")
        print('-' * 60)

        # this is teh end of the filename, the base filename will be added to it
        filenames_to_remove_list = ['_condition_builder_output.csv', '_final_result_dict_output.json',
                                    '_final_result_dict_output.pickle','_individual_log_and_result_files.zip',
                                    '_summary.txt']

        for partial_filename in filenames_to_remove_list:
            filename = f'{base_output_filename}{partial_filename}'

            if os.path.exists(filename):
                print(f"Deleting the file {filename}")
                os.remove(filename)

        # do we necessarily want to do this every time? removing folders will remove partially completed simulations...
        # now we need to remove any folders. we just remove any folder which has the basefilename_ in it, as taht will be a condition builder folder...

        remove_condition_builder_folders(mode='remove all folders', base_output_filename=base_output_filename, cwd='.')

    # now we need to be able to load the facility to get some information about it...
    facilities_folder = default_config_data['facilities_folder']

    facility_name = condition_builder_config_dict['facility']

    no_of_processes = condition_builder_config_dict['no_of_processes']

    print('-' * 60)
    print("Setting up the program.")
    print('-' * 60)

    print (f"Chosen facility is '{facility_name}'.")

    facility_yaml_filename = '{0}/{1}.yaml'.format(facilities_folder, facility_name)
    facility_yaml_file = open(os.path.expandvars(facility_yaml_filename))
    facility_input_data = yaml.load(facility_yaml_file , Loader=yaml.FullLoader)

    facility = Facility(facility_input_data)

    # get some values which we may need while working through the calculation
    facility_type = facility.get_facility_type()
    secondary_driver_flag = facility.get_secondary_driver_flag()
    nozzle_flag = facility.get_nozzle_flag()

    if 'driver_condition' in condition_builder_config_dict and 'driver_condition_list' in condition_builder_config_dict:
        raise Exception(f'pitot3_condition_builder() driver_condition and driver_condition_List cannot be in the same config file.')
    elif 'driver_condition' in condition_builder_config_dict:

        driver_condition_list = [condition_builder_config_dict['driver_condition']]

        base_config_dict.pop('driver_condition')

        print(f"Only the driver condition {condition_builder_config_dict['driver_condition']} will be tested.")

    elif 'driver_condition_list' in condition_builder_config_dict:

        driver_condition_list = condition_builder_config_dict['driver_condition_list']

        base_config_dict.pop('driver_condition_list')

        if len(driver_condition_list) == 1:
            print(f"Only the driver condition {driver_condition_list[0]} will be tested.")
        else:
            print(f"{len(driver_condition_list)} driver conditions will be tested: {driver_condition_list}")

    # we store the variables we are iterating through...
    variables_we_iterate_through = ['output_filename','driver_condition']

    # I put empty lists for the pressures here as then we can add to them below by calling them up in the namespace
    fill_states = []
    if secondary_driver_flag:
        fill_states += ['sd1']
        psd1_values = []

    fill_states += ['1']
    p1_values = []

    if facility_type == 'expansion_tube':
        fill_states += ['5']
        p5_values = []

    for fill_state in fill_states:

        if f'p{fill_state}' in condition_builder_config_dict and f'p{fill_state}_range' in condition_builder_config_dict:
            raise Exception(f'pitot3_condition_builder() p{fill_state} and p{fill_state}_range cannot be in the same config file.')
        elif f'p{fill_state}' in condition_builder_config_dict:

            fill_pressure = condition_builder_config_dict[f'p{fill_state}']

            locals()[f'p{fill_state}_values'] += [fill_pressure]

            base_config_dict.pop(f'p{fill_state}')

            print(f'p{fill_state} will be keep constant at {fill_pressure} Pa')

        elif f'p{fill_state}_range' in condition_builder_config_dict and f'p{fill_state}_no_of_samples' in condition_builder_config_dict \
                and f'p{fill_state}_spacing' in condition_builder_config_dict:

            fill_pressure_range = condition_builder_config_dict[f'p{fill_state}_range']
            fill_pressure_spacing = condition_builder_config_dict[f'p{fill_state}_spacing']
            fill_pressure_no_of_samples = condition_builder_config_dict[f'p{fill_state}_no_of_samples']

            if fill_pressure_spacing == 'linear':
                # trying to add our variable to the name space with a string title
                # put the variable in a dictionary and add the dictionary to the name space...
                locals()[f'p{fill_state}_values'] += list(np.linspace(fill_pressure_range[0], fill_pressure_range[1], fill_pressure_no_of_samples))

                fill_pressure_spacing_word = 'linearly'

            elif fill_pressure_spacing  == 'log':
                # trying to add our variable to the name space with a string title
                # put the variable in a dictionary and add the dictionary to the name space...
                locals()[f'p{fill_state}_values'] += list(np.geomspace(fill_pressure_range[0], fill_pressure_range[1], fill_pressure_no_of_samples))

                fill_pressure_spacing_word = 'logarithmically'

            base_config_dict.pop(f'p{fill_state}_range')
            base_config_dict.pop(f'p{fill_state}_spacing')
            base_config_dict.pop(f'p{fill_state}_no_of_samples')

            print(f'p{fill_state} will be varied from {fill_pressure_range[0]} to {fill_pressure_range[1]} Pa {fill_pressure_spacing_word} using {fill_pressure_no_of_samples} samples.')

        else:
            raise Exception(f'pitot3_condition_builder() Enough values have not been provided to make a p{fill_state} list. Need p{fill_state} or p{fill_state}_range, p{fill_state}_no_of_samples, AND p{fill_state}_spacing')

        variables_we_iterate_through += [f'p{fill_state}']

    # we multiply through by the other variables after...
    no_of_simulations_to_be_ran = len(driver_condition_list)
    if secondary_driver_flag:
        no_of_simulations_to_be_ran = no_of_simulations_to_be_ran*len(psd1_values)
    no_of_simulations_to_be_ran = no_of_simulations_to_be_ran*len(p1_values)
    if facility_type == 'expansion_tube':
        no_of_simulations_to_be_ran = no_of_simulations_to_be_ran*len(p5_values)
    print (f"Overall {no_of_simulations_to_be_ran} tests will be run.")

    #----------------------------------------------------------------------------------------
    # start by building the output folders...

    print('-'*60)
    print('Setting up the run folders for each simulation.')
    print('-'*60)

    test_names = []

    # we store the config that we are changing for each simulation in case we need it later on...
    changing_input_config_dict = {}

    lists_to_iterate_through = [driver_condition_list]

    if secondary_driver_flag:
        lists_to_iterate_through += [psd1_values]

    lists_to_iterate_through += [p1_values]

    if facility_type == 'expansion_tube':
        lists_to_iterate_through += [p5_values]

    for i, variables in enumerate(itertools.product(*lists_to_iterate_through)):

        driver_condition = variables[0]

        if secondary_driver_flag:
            psd1 = variables[1]
            p1 = variables[2]
            if facility_type == 'expansion_tube':
                p5 = variables[3]
        else:
            p1 = variables[1]
            if facility_type == 'expansion_tube':
                p5 = variables[2]

        test_number = i

        test_name = f'{base_output_filename}_{test_number}'

        test_names.append(test_name)

        run_folder = test_name

        pitot_3_input_file_filename = test_name

        config_dict = base_config_dict.copy()
        config_dict['output_filename'] = test_name
        config_dict['driver_condition'] = driver_condition
        if secondary_driver_flag:
            config_dict['psd1'] = psd1
        config_dict['p1'] = p1
        if facility_type == 'expansion_tube':
            config_dict['p5'] = p5
        config_dict['test_number'] = test_number

        pitot3_remote_run_creator(config_dict, run_folder, pitot_3_input_file_filename)

        # store the input variables in the dictionary we made for if we need it...

        changing_input_config_dict[test_name] = {}

        for variable in variables_we_iterate_through:
            changing_input_config_dict[test_name][variable] = config_dict[variable]

    # ----------------------------------------------------------------------------------------
    print('-' * 60)
    print('Checking for any excess old condition builder folders in the simulation folder and cleaning them up.')
    print('-' * 60)

    # this is important, as a previous simulation could have had more simulations...

    remove_condition_builder_folders(mode='conserve folders', base_output_filename=base_output_filename, folder_list=test_names, cwd='.')

    #----------------------------------------------------------------------------------------
    print('-'*60)
    print('Now starting to run the simulations.')
    print('-'*60)

    # get the start time here...

    start_time = time.perf_counter()

    # using some tips from here:
    # https://superfastpython.com/multiprocessing-pool-python/

    if no_of_processes == 1:
        for test_name in test_names:
            pitot3_condition_builder_test_run(test_name, changing_input_config_dict, variables_we_iterate_through)
    else:
        print(f"Automatically splitting calculations over {no_of_processes} cores.")
        print("This may make the terminal look a bit funny, so beware of this.")

        # multiprocessing!
        import multiprocessing

        pool = multiprocessing.Pool(processes = no_of_processes)

        # this will give us a list where only the test_name is changing...
        input_values = zip(test_names, itertools.repeat(changing_input_config_dict), itertools.repeat(variables_we_iterate_through))

        results = pool.starmap(pitot3_condition_builder_test_run, input_values)

        pool.close()

    end_time = time.perf_counter()

    total_time = end_time - start_time

    print('-'*60)
    print(f"This total calculation took {total_time:.2f} s ({total_time/60.0:.2f} mins, {total_time/3600.0:.2f} hours).")
    print("(Noting that this does not include the time taken for any runs performed before the last restart.)")
    print('-'*60)

    #----------------------------------------------------------------------------------------
    # when we're done, I guess we go through and make an output csv by re-loading the simulations and pulling out their one line outputs...
    # we will also build a results dictionary while we're at it...
    print('-'*60)
    print("Creating output csv and results dictionaries.")
    print('-'*60)

    # this will store a list for each variable from each test
    results_dict = {}

    # this will store the dict_of_objects for each simulation
    results_objects_dict = {}

    unsuccessful_simulations = []

    condition_builder_output_filename = f'{base_output_filename}_condition_builder_output.csv'

    with open(condition_builder_output_filename, 'w') as condition_builder_output_file:
        condition_builder_header = f"#Output of PITOT3 condition builder version {CONDITION_BUILDER_VERSION_STRING}"
        condition_builder_output_file.write(condition_builder_header + '\n')

        have_added_title_line = False

        for test_name in test_names:
            print('-'*60)
            print(f"Loading result from test {test_name}.")
            print('-'*60)

            # grab the test number too in case we need it ...
            test_number = test_name[len(base_output_filename)+1:]

            # where we start out...
            starting_working_directory = os.getcwd()

            # change directory to the one of the simulation
            os.chdir(starting_working_directory + '/' + test_name)

            files_in_the_current_run_directory = os.listdir(os.getcwd())

            json_filename = f'{test_name}.json'

            # this should ignore anything that failed...
            if  json_filename in files_in_the_current_run_directory:

                json_output_dict = pitot3_json_output_file_loader(json_filename)

                # I added the single line output stuff from each PITOT3 run to the .json output
                # so I could use it here instead of re-creating it like I did before...

                single_line_output_dict = json_output_dict['single_line_output_dict']

                title_line = single_line_output_dict['title_line']
                result_line = single_line_output_dict['result_line']
                title_list = single_line_output_dict['title_list']
                result_list = single_line_output_dict['result_list']

                if not have_added_title_line:
                    condition_builder_output_file.write(title_line + '\n')

                    have_added_title_line = True

                condition_builder_output_file.write(result_line + '\n')

                # add the result to the results dict
                for title, result in zip(title_list, result_list):
                    # add the variable if it isn't there yet...
                    if title not in results_dict:
                        results_dict[title] = []
                    results_dict[title].append(result)

                # add the json_output_dict to the dictionary for that too...

                results_objects_dict[test_name] = json_output_dict

            else:
                print(f"{test_name} does not have a .json output file so it must have failed.")
                unsuccessful_simulations.append(test_number)

            # return to the original directory when we're done... (this may actually be unnecessary except for at teh end?)
            os.chdir(starting_working_directory)

    #----------------------------------------------------------------------------------------
    # now we export the results dict to a json file...

    print('-'*60)
    print("Saving the simple results dictionary to a .json file.")
    print('-'*60)

    json_results_dict_output_filename = f'{base_output_filename}_final_result_dict_output.json'

    with open(json_results_dict_output_filename, "w") as output_file:
        json.dump(results_dict, output_file)

    #----------------------------------------------------------------------------------------
    # And we pickle the

    print('-'*60)
    print("Saving the json result for each simulation to a .pickle file.")
    print('-'*60)

    pickle_results_dict_output_filename = f'{base_output_filename}_final_result_dict_output.pickle'

    # may as well use the pitot3 pickle output function for this...
    pitot3_pickle_output_file_creator(results_objects_dict, pickle_results_dict_output_filename)

    #----------------------------------------------------------------------------------------
    # now zip up the result... this kind of works, but just needs some work as it also zips unnecessary stuff...

    print('-'*60)
    print("Zipping up the individual simulation results.")
    print('-'*60)

    cwd = '.'  # easier than the full path now...

    zipfile_name = f'{base_output_filename}_individual_log_and_result_files.zip'

    with zipfile.ZipFile(zipfile_name, "w") as zf:

        for dirname, subdirs, files in os.walk(cwd):

            #remove the working directory from the dirname and then check if what is left is in the test_names, if so,
            # we want to zip up that folder and all of its subfiles
            if len(dirname) > len(cwd): # we're not the working directory...
                dirname_for_comparison = dirname[len(cwd) + 1:] # +1 to remove the slash...

                if dirname_for_comparison in test_names:

                    zf.write(dirname)
                    for filename in files:
                        zf.write(os.path.join(dirname, filename))

    #----------------------------------------------------------------------------------------
    # now to finish off we delete the results folders that we just zipped up...

    print('-'*60)
    print("Now removing the individual simulation results to clean up the folder.")
    print('-'*60)

    remove_condition_builder_folders(mode = 'remove folders', folder_list = test_names, cwd = '.')

    #----------------------------------------------------------------------------------------
    # And make a summary of the simulation results and print it to the screen and to a file.

    print ('-'*60)
    print ("Printing summary of the total calculation to the screen and to a text document.")
    print ('-'*60)

    condition_builder_summary_filename = f"{base_output_filename}_condition_builder_summary.txt"

    with open(condition_builder_summary_filename, "w") as condition_builder_summary_file:

        output_list = [sys.stdout, condition_builder_summary_file]

        for output_stream in output_list:

            summary_line_1 = f'# Summary of PITOT3 condition builder. Performed using version {CONDITION_BUILDER_VERSION_STRING}.'

            print(summary_line_1, file=output_stream)

            no_of_successful_simulations = no_of_simulations_to_be_ran - len(unsuccessful_simulations)

            percentage_successful = (no_of_successful_simulations/no_of_simulations_to_be_ran)*100

            summary_line_2 = f"{no_of_simulations_to_be_ran} tests ran. {no_of_successful_simulations} ({percentage_successful :.1f}%) were successful."

            print(summary_line_2, file=output_stream)

            if unsuccessful_simulations:
                summary_line_3 = f"Unsuccessful simulations were test numbers {unsuccessful_simulations}"
                print(summary_line_3, file=output_stream)

            if facility_type == 'expansion_tube':
                facility_type_statement = 'an expansion tube'

            if nozzle_flag:
                nozzle_statement = 'with a nozzle'
            else:
                nozzle_statement = "without a nozzle"

            if secondary_driver_flag:
                secondary_driver_statement = ' and with a secondary driver.'
            else:
                secondary_driver_statement = '.'

            summary_line_4 = f'Calculations were performed for {facility_type_statement} {nozzle_statement}{secondary_driver_statement}'

            print(summary_line_4, file=output_stream)

            if len(driver_condition_list) > 1:
                summary_line_5 = f"{len(driver_condition_list)} driver conditions were tested ({driver_condition_list})."
            else:
                summary_line_5 = f"Only the driver condition {driver_condition_list[0]} was tested."

            print(summary_line_5, file=output_stream)

            variables_to_not_summarise_list = ['test_number', 'driver_condition','area_ratio',
                                               'secondary_driver_gas_gas_model', 'secondary_driver_gas_name',
                                               'test_gas_gas_model', 'test_gas_name',
                                               'accelerator_gas_gas_model', 'accelerator_gas_name']

            for variable in results_dict.keys():
                # now check the list and remove any string values from failed calculations
                data_list = []
                for value in results_dict[variable]:
                    if isinstance(value, (float, int)): data_list.append(value)
                if len(data_list) > 0:  # ie. don't bother summarising if there is no numerical data there
                    min_value = min(data_list)
                    max_value = max(data_list)

                    if variable[0] == 'p':
                        summary_line = f"Variable {variable} varies from {min_value:.2f} - {max_value:.2f} Pa."
                    elif variable == 'total_p':
                        summary_line = f"Variable {variable} varies from {min_value/1.0e6:.2f} - {max_value/1.0e6:.2f} MPa."
                    elif variable[0] in ['v','a'] or variable in ['Ue']:
                        summary_line = f"Variable {variable} varies from {min_value:.2f} - {max_value:.2f} m/s."
                    elif variable[0:2] == 'Ht' or variable[0] == 'h':
                        summary_line = f"Variable {variable} varies from {min_value/1.0e6:.2f} - {max_value/1.0e6:.2f} MJ/kg."
                    elif variable[0] == 'T':
                        summary_line = f"Variable {variable} varies from {min_value:.2f} - {max_value:.2f} K."
                    elif 'rho' in variable:
                        summary_line = f"Variable {variable} varies from {min_value:.4e} - {max_value:.4e} kg/m**3."
                    elif variable[0] == 'M' or 'gamma' in variable:
                        summary_line = f"Variable {variable} varies from {min_value:.2f} - {max_value:.2f}."
                    elif variable[0] == 'R':
                        summary_line = f"Variable {variable} varies from {min_value:.2f} - {max_value:.2f} J/kg.K."
                    elif 'unit_Re' in variable:
                        summary_line = f"Variable {variable} varies from {min_value:.4e} - {max_value:.4e} /m."
                    elif 'mu' in variable:
                        summary_line = f"Variable {variable} varies from {min_value:.4e} - {max_value:.4e} Pa.s."
                    elif variable == 'basic_test_time_us':
                        summary_line = f"Variable {variable} varies from {min_value:.2f} - {max_value:.2f} us."
                    elif 'X' in variable or 'c' in variable:
                        # this needs a bit more care as we get big and small numbers
                        if min_value == 0.0:
                            min_value_string = '0.0'
                        elif min_value < 0.01:
                            min_value_string = f'{min_value:.4e}'
                        else:
                            min_value_string = f'{min_value:.4f}'

                        if max_value == 0.0:
                            max_value_string = '0.0'
                        elif max_value < 0.01:
                            max_value_string = f'{max_value:.4e}'
                        else:
                            max_value_string = f'{max_value:.4f}'

                        summary_line = f"Variable {variable} varies from {min_value_string} - {max_value_string}."
                    else:
                        summary_line = f"Variable {variable} varies from {min_value:.3f} - {max_value:.3f}."

                    print(summary_line, file=output_stream)

    return

def pitot3_condition_builder_test_run(test_name, changing_input_config_dict, variables_we_iterate_through):
    """
    The function that runs the actual tests after the folders have been created.

    I made this a separate function to make it easy to parallelise it later on.

    :param test_name:
    :return:
    """

    print('-' * 60)
    print(f"Attempting to run test {test_name}.")
    print('-' * 60)

    # where we start out...
    starting_working_directory = os.getcwd()

    # change directory to the one of the simulation
    os.chdir(starting_working_directory + '/' + test_name)

    # now we want to look at what files we have in the current directory

    files_in_the_current_run_directory = os.listdir(os.getcwd())

    pitot_3_input_file_filename = f'{test_name}.yaml'

    # we implement a few checks here to see if we need to re-run the simulation or not

    if f'{test_name}.json' not in files_in_the_current_run_directory and f'{test_name}_failed.txt' not in files_in_the_current_run_directory:
        # this is one of the output files, so the simulation has not been ran.
        run_simulation = True
    elif f'{test_name}.json' not in files_in_the_current_run_directory and f'{test_name}_failed.txt' in files_in_the_current_run_directory:
        # {test_name}_failed.txt file is just a dummy file which tells the condition builder that this run failed.
        # it is no point re-running the simulation if this is in the folder!
        run_simulation = False
    else:
        # we need to load the json file to double check that this result is the same as the input file
        # (it could have been an old simulation)

        json_filename = f'{test_name}.json'

        json_output_dict = pitot3_json_output_file_loader(json_filename)

        config_data = json_output_dict['config_data']

        for variable in variables_we_iterate_through:
            if changing_input_config_dict[test_name][variable] != config_data[variable]:
                print(
                    f"Variable {variable} is different between the new simulation ({changing_input_config_dict[test_name][variable]}) and the saved one ({config_data[variable]}).")
                run_simulation = True
                break
            else:
                run_simulation = False

        if run_simulation:
            print("There was a result in this folder but at least one of the variables was different between it and the new simulation.")
            print("For this reason, a new simulation will be run.")

    if run_simulation:
        simulation_logfile = f'{test_name}.log'

        print(f"Piping run output to logfile {simulation_logfile}.")

        try:

            with open(simulation_logfile, 'w') as logfile:
                with redirect_stdout(logfile):
                    with redirect_stderr(logfile):
                        run_pitot3(config_filename=pitot_3_input_file_filename)
        except Exception as e:
            print("The run appears to have failed:")
            print(e)
            print("The result will not be added to the output.")

            failure_filename = f'{test_name}_failed.txt'
            print(f"Te dummy file {failure_filename} will also be created in the folder so the program knows that the run failed.")

            with open(failure_filename, 'w') as failure_file:
                pass
    else:
        print(f"Test {test_name} has already been ran. This test will not be re-ran.")

    # return to the original directory when we're done... (this may actually be uncessary?)
    os.chdir(starting_working_directory)

    return

def remove_condition_builder_folders(mode = 'remove all folders', base_output_filename = None, folder_list = None, cwd = '.'):
    """
    This is a function to remove condition builder folders.

    The user can provide the base filename and remove all folders with this name plus an underscore ('remove all folders' mode).

    They can also provide a folder list and CONSERVE these folders and delete any others ('conserve folders') or
    the user can provide a folder list and DELETE the specified folders ('delete folders')

    :return:
    """

    for dirname, subdirs, files in os.walk(cwd):

        #remove the working directory from the dirname and then check if what is left is in the test_names, if so,
        # we want to zip up that folder and all of its subfiles
        if len(dirname) > len(cwd): # we're not the working directory...
            dirname_for_comparison = dirname[len(cwd) + 1:] # +1 to remove the slash...

            if mode == 'remove all folders' and f'{base_output_filename}_' in dirname:
                    print(f"Deleting the folder {dirname}")
                    shutil.rmtree(dirname)

            elif mode == 'conserve folders' and dirname_for_comparison not in folder_list and f'{base_output_filename}_' in dirname:
                print(f"Deleting the folder {dirname} as it is not in the list of simulation folders to conserve.")
                shutil.rmtree(dirname)
            elif mode == 'remove folders' and dirname_for_comparison in folder_list:
                print(f"Deleting the folder {dirname} as it is in the list of folders to delete.")
                shutil.rmtree(dirname)

    return

# ----------------------------------------------------------------------------

def main():
    import optparse
    op = optparse.OptionParser(version=CONDITION_BUILDER_VERSION_STRING)
    op.add_option('-c', '--config_filename', '--config-filename', dest='config_filename',
                  help=("filename where the user configuration file is located"))

    opt, args = op.parse_args()
    config_filename = opt.config_filename

    run_pitot3_condition_builder(config_filename=config_filename)

# ----------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print ("PITOT3 condition builder  - Equilibrium impulse facility simulator")
        print ("start with --help for help with inputs")

    else:
        main()

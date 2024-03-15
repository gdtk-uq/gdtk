"""
pitot3_v2_loss_factor_optimiser.py:

This code allows the user to specify a measured Vs2 value and it will find the required V2 loss factor value to get that
Vs2 value in PITOT3. This run is fully theoretical, but a 'experimental_shock_tube_theoretical_acceleration_tube' version
(where vs1 is specified) could be made easily by changing a few things below,

Chris James (c.james4@uq.edu.au) - 15/03/24

"""

import numpy as np
import matplotlib.pyplot as plt

import time
import copy

import sys, os

from pitot3 import run_pitot3
from pitot3_utils.pitot3_classes import eilmer4_CEAGas_input_file_creator
from gdtk.gas import GasModel, GasState, GasFlow

from scipy.optimize import newton

start_time = time.perf_counter()

# This is kind of messing around with the 12A condition from Curtis et al. (2022)
# Optimisation of super orbital Earth re-entry conditions in an expansion tube

config_dict = {'mode': 'fully_theoretical', 'output_filename': f'loss_factor_test',
               'facility': 'x2_nozzle', 'driver_condition': 'x2-lwp-2.5mm-100He-0-empirical',
               'test_gas_gas_model': 'CEAGas', 'test_gas_name': 'n2-o2-with-ions', 'p1': 1800.0, 'p5':13.2,
               'secondary_diaphragm_diaphragm_type' : 'velocity_loss_factor', 'secondary_diaphragm_velocity_loss_factor':1.0}

vs2_experimental = 12000.0 # m/s

print('-' * 60)
print(f"Experimental vs2 = {vs2_experimental} m/s.")

def error_in_vs2_eqn(velocity_loss_factor, pitot3_config_dict = config_dict, vs2_experimental = vs2_experimental):
    print('-' * 60)

    print(f"Guessed velocity loss factor = {velocity_loss_factor:.6f}.")

    print("Running pitot 3 with this velocity loss factor.")

    # doing a deepcopy here so we can change values but not mess with the original values...
    pitot3_config_dict = copy.deepcopy(pitot3_config_dict)

    # change the compression ratio to the guessed one...
    pitot3_config_dict['secondary_diaphragm_velocity_loss_factor'] = velocity_loss_factor

    # and the output filename too...

    # we make the original filename a folder... (it will only be created the first time)

    os.makedirs(pitot3_config_dict['output_filename'], exist_ok=True)

    pitot3_config_dict['output_filename'] += f'/velocity_loss_factor_{velocity_loss_factor:.2f}'

    config_data, gas_path, object_dict, states_dict = run_pitot3(config_dict = pitot3_config_dict)

    acceleration_tube_object = object_dict['acceleration_tube']

    vs2 = acceleration_tube_object.get_shock_speed()

    print('-' * 60)

    print(f"Found vs2 = {vs2} m/s.")

    print(f"vs2 found - vs2 experimental = {vs2 - vs2_experimental} m/s.")

    return vs2 - vs2_experimental

# just use the ideal value from above...
velocity_loss_factor_guess = config_dict['secondary_diaphragm_velocity_loss_factor']

# second guess here is just to make sure it goes smaller
velocity_loss_factor, results_info = newton(error_in_vs2_eqn, x0 = velocity_loss_factor_guess, x1=velocity_loss_factor_guess*0.99,
                                            tol=1.0e-3, full_output=True)

iterations = results_info.iterations

end_time = time.perf_counter()

calculation_time = end_time - start_time

text_output_filename = config_dict['output_filename'] + '_final_text_output.txt'

with open(text_output_filename, "w") as output_file:

    for output_stream in [sys.stdout, output_file]:

        print('-' * 60, file=output_stream)
        print(f"Experimental vs2 was {vs2_experimental} m/s.", file=output_stream)
        print(f"Final found velocity loss factor is {velocity_loss_factor:.6f}.", file=output_stream)
        print(f"Total calculation took {calculation_time:.2f} s ({calculation_time/60.0:.2f} minutes) and {iterations} iterations.", file=output_stream)
        print('-' * 60, file=output_stream)



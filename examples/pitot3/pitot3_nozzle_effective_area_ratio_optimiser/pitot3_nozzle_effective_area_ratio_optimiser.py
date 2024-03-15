"""
pitot3_nozzle_effective_area_ratio_optimiser.py

Given a PITOT3 config dictionary, this code will find the effective area ratio which reaches the desired Pitot pressure.
This is a theoretical case, which I doubt one would use, but it could be used with an experimental simulation to
get the effective nozzle area ratio for a measured Pitot pressure.

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

output_filename = 'nozzle_effective_area_ratio_results'

# This is kind of messing around with the 12A condition from Curtis et al. (2022)
# Optimisation of super orbital Earth re-entry conditions in an expansion tube

config_dict = {'mode': 'fully_theoretical', 'output_filename': output_filename,
               'facility': 'x2_nozzle', 'driver_condition':'x2-lwp-2.5mm-100He-0-empirical',
               'test_gas_gas_model': 'CEAGas', 'test_gas_name': 'n2-o2-with-ions', 'p1': 1800.0, 'p5':13.2,
               'secondary_diaphragm_diaphragm_type' : 'velocity_loss_factor', 'secondary_diaphragm_velocity_loss_factor':1.0}

desired_pitot_pressure = 200.0e3 # Pa

area_ratio_guess_1 = 5.64 # geometric
area_ratio_guess_2 = 5.0 #and a bit less

print('-' * 60)
print(f"Desired pitot_pressure = {desired_pitot_pressure/1.0e3:.2f} kPa ({desired_pitot_pressure:.2f} Pa).")

def error_in_pitot_pressure_eqn(area_ratio, pitot3_config_dict = config_dict, desired_pitot_pressure = desired_pitot_pressure):
    print('-' * 60)

    print(f"Guessed area_ratio  = {area_ratio:.2f}")

    print("Running PITOT3 with this area_ratio value.")

    # doing a deepcopy here so we can change values but not mess with the original values...
    pitot3_config_dict = copy.deepcopy(pitot3_config_dict)

    # change the compression ratio to the guessed one...
    pitot3_config_dict['area_ratio'] = area_ratio

    # and the output filename too...

    # we make the original filename a folder... (it will only be created the first time)

    os.makedirs(pitot3_config_dict['output_filename'], exist_ok=True)

    pitot3_config_dict['output_filename'] += f'/effective_area_ratio_{area_ratio:.2f}'

    config_data, gas_path, object_dict, states_dict = run_pitot3(config_dict = pitot3_config_dict)

    pitot_pressure = states_dict['s8'].get_pitot_pressure()

    print('-' * 60)

    print(f"Found pitot_pressure = {pitot_pressure/1.0e3:.2f} kPa ({pitot_pressure:.2f}) Pa.")

    print(f"pitot_pressure found - desired pitot_pressure = {(pitot_pressure - desired_pitot_pressure)/1.0e3:.4f} kPa ({(pitot_pressure - desired_pitot_pressure):.4f} Pa).")

    return pitot_pressure - desired_pitot_pressure

# second guess here is just to make sure it goes smaller
area_ratio, results_info = newton(error_in_pitot_pressure_eqn, x0 = area_ratio_guess_1, x1=area_ratio_guess_2,
                                            tol=1.0e-2, full_output=True)

iterations = results_info.iterations

end_time = time.perf_counter()

calculation_time = end_time - start_time

text_output_filename = config_dict['output_filename'] + '_final_text_output.txt'

with open(text_output_filename, "w") as output_file:

    for output_stream in [sys.stdout, output_file]:

        print('-' * 60, file=output_stream)
        print(f"Desired pitot_pressure was {desired_pitot_pressure/1.0e3:.2f} kPa ({desired_pitot_pressure:.2f} Pa).", file=output_stream)
        print(f"Final found area_ratio is {area_ratio:.2f}.", file=output_stream)
        print(f"Total calculation took {calculation_time:.2f} s ({calculation_time/60.0:.2f} minutes) and {iterations} iterations.", file=output_stream)
        print('-' * 60, file=output_stream)



"""
pitot3_x2_effective_driver_condition_optimiser.py

This example uses PITOT3 with a basic optimiser to allow us to work out the experimental compression ratio
from a measured experimental Vs1 value. This is a bit of an assumption as it assumes that the main variable controlling
driver performance is the compression ratio, but I think it is a pretty good ones, as that controls the temperature,
which controls the driver sound speed, and with it, the performance.

It guesses using the original compression ratio which should be pretty good.

It runs in shock tube mode as the driver only affects Vs1. This is obviously fine for simulating the X2 driver.

There is a similar, slightly more complicated procedure in Gildfind et al. (2015)
Free-piston driver performance characterisation using experimental shock speeds through helium
for those who are interested.

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

# just noting that this is not a real condition, I was just playing around...

# default driver values, we will change compression ratio to match stuff up...
driver_dict = {'driver_condition_name': 'x2-lwp-2.5mm-0', 'driver_condition_type': 'isentropic-compression-compression-ratio',
               'driver_gas_model': 'CEAGas',
               'driver_fill_composition': {'He': 0.8, 'Ar': 0.2}, 'driver_speciesList': ['He', 'Ar'],
               'driver_inputUnits': 'moles', 'driver_withIons': False,
               'driver_p': 77.2e3, 'driver_T': 298.15, 'compression_ratio': 40.0,
               'D_throat': 0.085}

# just nrst mode so it runs fast
config_dict = {'mode': 'fully_theoretical', 'output_filename': 'compression_ratio_test_results',
               'facility': 'x2_nrst_85_mm_shock_tube', 'driver_condition': 'custom_from_dict',
               'driver_dict': driver_dict,
               'test_gas_gas_model': 'CEAGas', 'test_gas_name': 'n2-o2-with-ions', 'p1': 1800.0}

vs1_experimental = 5200.0 # m/s

print('-' * 60)
print(f"Experimental vs1 = {vs1_experimental} m/s.")

def error_in_vs1_eqn(compression_ratio, pitot3_config_dict = config_dict, vs1_experimental = vs1_experimental):
    print('-' * 60)

    print(f"Guessed compression ratio = {compression_ratio:.2f}.")

    print("Running pitot 3 with this compression ratio.")

    # doing a deepcopy here so we can change values but not mess with the original values...
    pitot3_config_dict = copy.deepcopy(pitot3_config_dict)

    # change the compression ratio to the guessed one...
    pitot3_config_dict['driver_dict']['compression_ratio'] = compression_ratio

    # and the output filename too...

    # we make the original filename a folder... (it will only be created the first time)

    os.makedirs(pitot3_config_dict['output_filename'], exist_ok=True)

    pitot3_config_dict['output_filename'] += f'/compression_ratio_{compression_ratio:.2f}'

    config_data, gas_path, object_dict, states_dict = run_pitot3(config_dict = pitot3_config_dict)

    shock_tube_object = object_dict['shock_tube']

    vs1 = shock_tube_object.get_shock_speed()

    print('-' * 60)

    print(f"Found vs1 = {vs1} m/s.")

    print(f"vs1 found - vs1 experimental = {vs1 - vs1_experimental} m/s.")

    return vs1 - vs1_experimental

# just use the ideal value from above...
compression_ratio_guess = driver_dict['compression_ratio']

compression_ratio, results_info = newton(error_in_vs1_eqn, compression_ratio_guess, tol=1.0e-2, full_output=True)

iterations = results_info.iterations

end_time = time.perf_counter()

calculation_time = end_time - start_time

text_output_filename = config_dict['output_filename'] + '_final_text_output.txt'

with open(text_output_filename, "w") as output_file:

    for output_stream in [sys.stdout, output_file]:
        print('-' * 60, file=output_stream)
        print(f"Specified Experimental vs1 was {vs1_experimental} m/s.", file=output_stream)
        print(f"Final found compression ratio is {compression_ratio:.2f}.", file=output_stream)
        print(f"Total calculation took {calculation_time:.2f} s ({calculation_time/60.0:.2f} minutes) and {iterations} iterations.", file=output_stream)
        print('-' * 60, file=output_stream)



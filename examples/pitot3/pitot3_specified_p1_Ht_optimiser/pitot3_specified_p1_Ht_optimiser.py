"""
pitot3_specified_p1_Ht_optimiser.py:

For a given shock tube case (driver, p1, test gas etc.) this code will find the p5 that gives the required enthalpy.

If you give it bad first guesses for p5 it will take a while to solve, if the guesses are good, it'll be very quick,
but there is an element of chance there!

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

config_dict = {'mode': 'fully_theoretical', 'output_filename': 'Ht_optimiser_results',
               'facility': 'x2_nozzle', 'driver_condition':'x2-lwp-2.5mm-100He-0-empirical',
               'test_gas_gas_model': 'CEAGas', 'test_gas_name': 'n2-o2-with-ions', 'p1': 1800.0,
               'secondary_diaphragm_diaphragm_type' : 'velocity_loss_factor', 'secondary_diaphragm_velocity_loss_factor':1.0}

desired_Ht = 85.0e6 # J/kg

p5_guess_1 = 10.0 # Pa
p5_guess_2 = 9.0 # Pa

print('-' * 60)
print(f"Desired Ht = {desired_Ht/1.0e6:.2f} MJ/kg ({desired_Ht:.2e} J/kg).")

def error_in_Ht_eqn(p5, pitot3_config_dict = config_dict, desired_Ht = desired_Ht):
    print('-' * 60)

    print(f"Guessed p5 = {p5:.6f} Pa.")

    print("Running PITOT3 with this p5 value.")

    # doing a deepcopy here so we can change values but not mess with the original values...
    pitot3_config_dict = copy.deepcopy(pitot3_config_dict)

    # change the compression ratio to the guessed one...
    pitot3_config_dict['p5'] = p5

    # and the output filename too...

    # we make the original filename a folder... (it will only be created the first time)

    os.makedirs(pitot3_config_dict['output_filename'], exist_ok=True)

    pitot3_config_dict['output_filename'] += f'/p5_{p5}_Pa'

    config_data, gas_path, object_dict, states_dict = run_pitot3(config_dict = pitot3_config_dict)

    Ht = states_dict['s8'].get_total_enthalpy()

    print('-' * 60)

    print(f"Found Ht = {Ht/1.0e6:.2f} MJ/kg ({Ht:.2e} J/kg).")

    print(f"Ht found - desired Ht = {(Ht - desired_Ht)/1.0e6} MJ/kg ({Ht - desired_Ht:.4e} J/kg).")

    return Ht - desired_Ht


# second guess here is just to make sure it goes smaller
p5, results_info = newton(error_in_Ht_eqn, x0 = p5_guess_1, x1=p5_guess_2,
                                            tol=1.0e-2, full_output=True)

iterations = results_info.iterations

end_time = time.perf_counter()

calculation_time = end_time - start_time

text_output_filename = config_dict['output_filename'] + '_final_text_output.txt'

with open(text_output_filename, "w") as output_file:

    for output_stream in [sys.stdout, output_file]:

        print('-' * 60, file=output_stream)
        print(f"Desired Ht was {desired_Ht/1.0e6:.2f} MJ/kg ({desired_Ht:.2e} J/kg).", file=output_stream)
        print(f"Final found p5 is {p5:.6f} Pa.", file=output_stream)
        print(f"Total calculation took {calculation_time:.2f} s ({calculation_time/60.0:.2f} minutes) and {iterations} iterations.", file=output_stream)
        print('-' * 60, file=output_stream)



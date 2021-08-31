"""
x2_zander_condition_air_theory_scripting_example.py

This is the same as the x2_zander_condition_air_theory
example except now the simulation is specified as a
Python script which can be ran by going:

$ python3 x2_zander_condition_air_theory_scripting_example.py

It is hoped that this example will serve as a springboard for
the automation of PITOT3 in various programs.

Chris James (c.james4@uq.edu.au) - 30/08/21

"""

from pitot3 import run_pitot3

config_dict = {'mode':'fully_theoretical','output_filename':'x2_zander_condition_air_theory',
               'facility':'x2_nozzle', 'driver_condition':'x2-lwp-2.0mm-0',
               'test_gas_gas_model':'CEAGas', 'test_gas_name':'n2-o2-with-ions', 'p1':3000.0, 'p5':'10.0',
               'area_ratio':5.64, 'cone_half_angle_degrees':15.0, 'wedge_angle_degrees':54.0}

config_data, gas_path, object_dict = run_pitot3(config_dict = config_dict)







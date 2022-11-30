"""
x2_zander_condition_air_theory_scripting_example.py

This is the same as the x2_zander_condition_air_theory
example except now the simulation is specified as a
Python script which can be ran by going:

$ python3 x2_zander_condition_air_theory_scripting_example.py

It is hoped that this example will serve as a springboard for
the automation of PITOT3 in various programs.

I recently added some extra code at the bottom to show how to pull out some common values people might want - 29/11/22

Chris James (c.james4@uq.edu.au) - 30/08/21

"""

#----------------------------------------------------------------------------------------------
# the code that runs PITOT3

from pitot3 import run_pitot3

config_dict = {'mode':'fully_theoretical','output_filename':'x2_zander_condition_air_theory',
               'facility':'x2_nozzle', 'driver_condition':'x2-lwp-2.0mm-0',
               'test_gas_gas_model':'CEAGas', 'test_gas_name':'n2-o2-with-ions', 'p1':3000.0, 'p5':'10.0',
               'area_ratio':5.64, 'cone_half_angle_degrees':15.0, 'wedge_angle_degrees':54.0}

config_data, gas_path, object_dict, states_dict = run_pitot3(config_dict = config_dict)

#----------------------------------------------------------------------------------------------
# the code that lets us pull some useful things out of the output of the code

# the object dictionary has all of the states we need
# 'driver' is the driver
# 'shock_tube' is shock tube
# 'acceleration_tube' is the acceleration tube
# 'nozzle' is the nozzle
# 'test_section' is the test section
# there are also states for the diaphragms etc. as well

driver_object = object_dict['driver']
shock_tube_object = object_dict['shock_tube']
acceleration_tube_object = object_dict['acceleration_tube']
test_section_object = object_dict['test_section']

# lets us start by getting the shock speeds:

vs1 = shock_tube_object.get_shock_speed()
vs2 = acceleration_tube_object.get_shock_speed()

print('-'*60)
print("The calculated shock speeds are:")
print(f"vs1 = {vs1:.2f} m/s")
print(f"vs2 = {vs2:.2f} m/s")

# it is useful to show how to get a tube fill state, shocked state and unsteadily expanded state
# so I'll get the shock tube fill pressure (state1), shock tube post-shock state (state2),
# and the unsteadily expanded state in the acceleration tube (state 7)

state1_facility_state = shock_tube_object.get_fill_state()
state1_gas_state = state1_facility_state.get_gas_state()
state1_fill_pressure = state1_gas_state.p
state1_fill_temperature = state1_gas_state.T

print('-'*60)
print(f"Shock tube fill pressure (p1) is {state1_fill_pressure:.2f} Pa")
print(f"Shock tube fill temperature (T1) is {state1_fill_temperature:.2f} K")

state2_facility_state = shock_tube_object.get_shocked_state()
state2_gas_state = state2_facility_state.get_gas_state()
state2_pressure = state2_gas_state.p
state2_temperature = state2_gas_state.T
state2_composition, state2_composition_units = state2_facility_state.get_reduced_composition_with_units()

print('-'*60)
print(f"Shock tube post-shock pressure (p2) is {state2_pressure:.2f} Pa")
print(f"Shock tube post-shock temperature (T2) is {state2_temperature:.2f} K")
print(f"Shock tube post-shock composition (by {state2_composition_units}):")
print(state2_composition)

state7_facility_state = acceleration_tube_object.get_unsteadily_expanded_state()
state7_gas_state = state7_facility_state.get_gas_state()
state7_velocity = state7_facility_state.get_v()
state7_Mach_number = state7_facility_state.get_M()

print('-'*60)
print(f"Acceleration tube unsteadily expanded gas velocity (v7) is {state7_velocity:.2f} m/s.")
print(f"Acceleration tube unsteadily expanded gas Mach number (M7) is {state7_Mach_number:.2f}.")

# and then it is useful to get some test section information
# so we'll get some information about the test section freestream state (state8), and the test section post-shock
# frozen and equilibrium states(state 10f and state 10e)

state8_facility_state = test_section_object.get_test_section_state()
state8_gas_state = state8_facility_state.get_gas_state()
state8_pressure = state8_gas_state.p
state8_temperature = state8_gas_state.T
state8_density = state8_gas_state.rho
state8_gamma = state8_gas_state.gamma
state8_R = state8_gas_state.R
state8_velocity = state8_facility_state.get_v()
state8_Mach_number = state8_facility_state.get_M()

state8_composition, state8_composition_units = state8_facility_state.get_reduced_composition_with_units()

# add some extra values we might want, Ht, Ue, pitot pressure etc.

Ht = state8_facility_state.get_total_enthalpy()
Ue = state8_facility_state.get_flight_equivalent_velocity()
pitot_pressure = state8_facility_state.get_pitot_pressure()

print('-'*60)
print(f"Test section freestream pressure (p8) is {state8_pressure:.2f} Pa")
print(f"Test section freestream temperature (T8) is {state8_temperature:.2f} K")
print(f"Test section freestream density (rho8) is {state8_density:.2e} kg/m**3")
print(f"Test section freestream gamma and R are {state8_gamma:.2f} and {state8_R:.2f} J/kg/K")
print(f"Test section freestream velocity (v8) is {state8_velocity:.2f} m/s.")
print(f"Test section freestream Mach number (M8) is {state8_Mach_number:.2f}.")
print(f"Test section freestream composition (by {state8_composition_units}):")
print(state8_composition)

print(f"Test section freestream stagnation enthalpy (Ht) is {Ht/1.0e6:.2f} MJ/kg")
print(f"Test section freestream flight equivalent velocity (Ht) is {Ue:.2f} m/s")
print(f"Test section freestream Pitot pressure is {pitot_pressure/1000.0:.2f} kPa")

state10f_facility_state = test_section_object.get_post_normal_shock_ideal_gas_state()
state10f_gas_state = state10f_facility_state.get_gas_state()
state10f_pressure = state10f_gas_state.p
state10f_temperature = state10f_gas_state.T

print('-'*60)
print(f"Test section frozen post-normal shock pressure (p10f) is {state10f_pressure:.2f} Pa")
print(f"Test section frozen post-normal shock temperature (T10f) is {state10f_temperature:.2f} K")

state10e_facility_state = test_section_object.get_post_normal_shock_state()
state10e_gas_state = state10e_facility_state.get_gas_state()
state10e_pressure = state10e_gas_state.p
state10e_temperature = state10e_gas_state.T
state10e_composition, state10e_composition_units = state10e_facility_state.get_reduced_composition_with_units()

print('-'*60)
print(f"Test section equilibrium post-normal shock pressure (p10e) is {state10e_pressure:.2f} Pa")
print(f"Test section equilibrium post-normal shock temperature (T10e) is {state10e_temperature:.2f} K")
print(f"Test section equilibrium post-normal shock composition (by {state10e_composition_units}):")
print(state10e_composition)

















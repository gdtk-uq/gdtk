"""
pitot3_classes.py

Just a file to store all of the classes for PITOT3,
to make the main file a bit cleaner.

Chris James (c.james4@uq.edu.au) 01/01/21

"""

import sys, os, math
import yaml

from eilmer.gas import GasModel, GasState, GasFlow
from eilmer.ideal_gas_flow import p0_p
from eilmer.zero_solvers import secant

def eilmer4_CEAGas_input_file_creator(output_filename, mixtureName, speciesList, reactants,
                                      inputUnits, withIons, trace = 1.0e-6):
    """Just a function to make an input file for the CEAGas object..."""

    with open('{0}.lua'.format(output_filename), 'w') as gas_file:
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

    def __init__(self, cfg, p_0 = None, T_0 = None):
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

        self.M_throat = cfg['M_throat']

        if self.M_throat > 0.0:
            self.driver_exit_state_name = '3s'
        else:
            self.driver_exit_state_name = '4'

        # -------------------- gas composition ----------------------------

        self.driver_gas_model = cfg['driver_gas_model']

        # just doing the CEA model for now...
        if self.driver_gas_model == 'CEAGas':

            self.driver_fill_composition = cfg['driver_fill_composition']
            self.driver_speciesList = cfg['driver_speciesList']
            self.driver_inputUnits = cfg['driver_inputUnits']
            self.driver_withIons = cfg['driver_withIons']

            # now build a gas model file and attach it to this object...

            eilmer4_CEAGas_input_file_creator('PITOT3_cea_driver_condition', 'driver_gas', self.driver_speciesList, self.driver_fill_composition,
                                              self.driver_inputUnits, self.driver_withIons)

            self.gmodel = GasModel('PITOT3_cea_driver_condition.lua')

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
                print("Setting driver temperature to the default value of {0:.2f} K".format(T))
                self.driver_T = T_0

            # we call this state 4i
            state4i = GasState(self.gmodel)
            state4i.p = self.driver_p
            state4i.T = self.driver_T

            state4i.update_thermo_from_pT()
            state4i.update_sound_speed()

            gamma = state4i.gamma

            # we assume that the state 4 driver is stationary
            v4i = 0.0

            # make a reference gas state here...
            reference_gas_state = GasState(self.gmodel)
            reference_gas_state.p = p_0
            reference_gas_state.T = T_0

            reference_gas_state.update_thermo_from_pT()
            reference_gas_state.update_sound_speed()

            # now make our facility driver object...
            self.state4i = Facility_State('s4i', state4i, v4i, reference_gas_state=reference_gas_state)

            if self.driver_condition_type == 'isentropic-compression-p4':
                self.p4 = float(cfg['p4'])

                print ("Performing isentropic compression from the driver fill condition to {0:.2f} MPa.".format(self.p4/1.0e6))

                self.T4 = state4i.T * (self.p4 / state4i.p) ** (1.0 - (1.0 / gamma))  # K

            elif self.driver_condition_type == 'isentropic-compression-compression-ratio':

                self.compression_ratio = cfg['compression_ratio']

                print ("Performing isentropic compression from driver fill condition over compression ratio of {0}.".format(cfg['compression_ratio']))
                pressure_ratio =  self.compression_ratio**gamma #pressure ratio is compression ratio to the power of gamma

                self.p4 = state4i.p*pressure_ratio

            self.T4 = state4i.T * (self.p4 / state4i.p) ** (1.0 - (1.0 / gamma))  # K

        state4 = GasState(self.gmodel)
        state4.p = self.p4
        state4.T = self.T4

        state4.update_thermo_from_pT()
        state4.update_sound_speed()

        # we assume that the state 4 driver is stationary
        v4 = 0.0

        # make a reference gas state here...
        reference_gas_state = GasState(self.gmodel)
        reference_gas_state.p = p_0
        reference_gas_state.T = T_0

        reference_gas_state.update_thermo_from_pT()
        reference_gas_state.update_sound_speed()

        # now make our facility driver object...
        self.state4 = Facility_State('s4', state4, v4,
                                     reference_gas_state=reference_gas_state)

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
                                          reference_gas_state=reference_gas_state)

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
        print("Doing reflected shock at the {0} which the user has asked for.".format(self.diaphragm_name))

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
            print("Performing a reflected shock with a user selected Mach number of {0}".format(self.Mr_input))
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
                                                      reference_gas_state=reference_gas_state)

        print('Mr = {0:.2f}, vr = {1:.2f}'.format(self.Mr, self.vr))
        print(self.reflected_shocked_state)
        if self.reflected_shocked_state.get_gas_state().gmodel.type_str == 'CEAGas':
            print('species in {0} at equilibrium (by mass):'.format(self.reflected_shocked_state.get_state_name()))
            print(self.reflected_shocked_state.get_reduced_species_massf_dict())

        return

    def perform_velocity_loss_factor_at_diaphragm(self):

        print('-' * 60)
        print("Multiplying the velocity at the {0} by a velocity_loss_factor of {1}".format(self.diaphragm_name,
                                                                                            self.velocity_loss_factor))

        entrance_gas_gas_state = self.diaphragm_entrance_state.get_gas_state()
        entrance_gas_v = self.diaphragm_entrance_state.get_v()

        if self.diaphragm_entrance_state.reference_gas_state:
            reference_gas_state = self.diaphragm_entrance_state.get_reference_gas_state()
        else:
            reference_gas_state = None

        self.velocity_loss_factor_state = Facility_State(self.diaphragm_entrance_state_name + 'l',
                                                         entrance_gas_gas_state,
                                                         entrance_gas_v * self.velocity_loss_factor,
                                                         reference_gas_state=reference_gas_state)
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

    def __init__(self, state_name, gas_state, v, reference_gas_state = None):
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

        return

    def __str__(self):
        #text = "Facility_State(state {0}: p = {1:.2f} Pa, T = {2:.2f} K, V = {3:.2f} m/s, M = {4:.2f})".format(self.state_name,
        #                                                                                                       self.gas_state.p,
        #                                                                                                       self.gas_state.T,
        #                                                                                                       self.v, self.M)
        text = "Facility_State(state {0}: p = {1:.2f} Pa, T = {2:.2f} K, gam = {3:.2f}, R = {4:.2f} J/kg K, V = {5:.2f} m/s, M = {6:.2f})".format(self.state_name,
                                                                                                               self.gas_state.p,
                                                                                                               self.gas_state.T,
                                                                                                               self.gas_state.gamma,
                                                                                                               self.gas_state.R,
                                                                                                               self.v, self.M)

        return text

    def calculate_pitot_and_total_conditions(self):
        """
        Function to calculate the pitot and total conditions for the state, which can then be returned
        using the related functions.

        Remember that these states don't need the velocity object as they have been brought to rest.
        :return:
        """

        # start by making a gas flow object which we can use...
        gas_flow_object = GasFlow(self.gas_state.gmodel)

        # TODO: add stuff to check that these actually worked, and not sure them if they didn't...

        pitot_state = GasState(self.gas_state.gmodel)
        gas_flow_object.pitot_condition(self.gas_state, self.v, pitot_state)
        self.pitot_state = pitot_state

        total_state = GasState(self.gas_state.gmodel)
        gas_flow_object.total_condition(self.gas_state, self.v, total_state)
        self.total_state = total_state

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

    def get_pitot_condition(self):
        """
        Returns the pitot condition. Will calculate it if it doesn't exist yet.
        :return:
        """

        if hasattr(self, 'pitot_state'):
            return self.pitot_state
        else:
            #print("Will calculate the pitot condition and then return it.")
            self.calculate_pitot_and_total_conditions()
            return self.pitot_state

    def get_pitot_pressure(self):
        """
        Returns the pitot pressure (the pressure for the pitot condition). Will calculate it if it doesn't exist yet.
        :return:
        """

        if hasattr(self, 'pitot_state'):
            return self.pitot_state.p
        else:
            #print("Will calculate the pitot condition and then return it.")
            self.calculate_pitot_and_total_conditions()
            return self.pitot_state.p

    def get_total_condition(self):
        """
        Returns the total condition. Will calculate it if it doesn't exist yet.
        :return:
        """

        if hasattr(self, 'total_state'):
            return self.total_state
        else:
            self.calculate_pitot_and_total_conditions()
            return self.total_state

    def get_total_pressure(self):
        """
        Returns the total pressure (the pressure of the total condition). Will calculate it if it doesn't exist yet.
        :return:
        """

        if hasattr(self, 'total_state'):
            return self.total_state.p
        else:
            self.calculate_pitot_and_total_conditions()
            return self.total_state.p

    def get_total_temperature(self):
        """
        Returns the total temperature (the temperature of the total condition). Will calculate it if it doesn't exist yet.
        :return:
        """

        if hasattr(self, 'total_state'):
            return self.total_state.T
        else:
            self.calculate_pitot_and_total_conditions()
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

            reference_gas_state = self.get_reference_gas_state()

            # total enthalpy is the sensible enthalpy of the total condition - the sensible enthalpy of the reference state
            total_enthalpy = total_state.enthalpy - reference_gas_state.enthalpy

            self.total_enthalpy = total_enthalpy

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

        self.flight_equivalent_velocity = math.sqrt(2.0 * self.total_enthalpy)

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

    def get_reduced_species_massf_dict(self):
        """
        If the gas is a CEAGas, it will go through the gas dictionary and remove any values which are empty,
        which is very easy for printing things like fill states which probably only have 1-3 gases out of the full set of
        high temperature possibilities.

        Obviously not useful if the GasModel isn't a CEAGas...

        """

        if self.get_gas_state().gmodel.type_str == 'CEAGas':
            # a fill state will not have many different species, so we should just the species which actually exist
            species_massf_dict = self.get_gas_state().ceaSavedData['massf']
            reduced_species_massf_dict = {}
            for species in species_massf_dict.keys():
                if species_massf_dict[species] > 0.0:
                    reduced_species_massf_dict[species] = species_massf_dict[species]

            return reduced_species_massf_dict

        else:
            print("GasModel is not a CEAGas, so this function isn't useful. Will return None.")
            return None

    # def get_reduced_species_molef_dict(self):
    #     """
    #     If the gas is a CEAGas, it will go through the gas dictionary and remove any values which are empty,
    #     which is very easy for printing things like fill states which probably only have 1-3 gases out of the full set of
    #     high temperature possibilities.
    #
    #     Obviously not useful if the GasModel isn't a CEAGas...
    #
    #     This is the same as the function above but it converts to mole fractions while doing it...
    #
    #     """
    #
    #     if self.get_gas_state().gmodel.type_str == 'CEAGas':
    #         # a fill state will not have many different species, so we should just the species which actually exist
    #         species_massf_dict = self.get_gas_state().ceaSavedData['massf']
    #
    #         reduced_species_massf_dict = {}
    #         for species in species_massf_dict.keys():
    #             if species_massf_dict[species] > 0.0:
    #                 reduced_species_massf_dict[species] = species_massf_dict[species]
    #
    #         return reduced_species_massf_dict
    #
    #     else:
    #         print("GasModel is not a CEAGas, so this function isn't useful. Will return None.")
    #         return None

    def get_mu(self):
        """
        Function to return the dynamic viscosity mu, as it seems that one needs to the update the transport
        coefficients for it to exist at all...

        Mu is in Pa.s
        :return:
        """

        gas_state = self.get_gas_state()
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
                 preset_gas_models_folder, unsteady_expansion_steps, vs_guess_1, vs_guess_2, vs_limits, vs_tolerance):

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
        elif self.fill_gas_model == 'custom' and self.fill_gas_filename:
            fill_gmodel_location = self.fill_gas_filename

        fill_gmodel = GasModel(os.path.expandvars(fill_gmodel_location))

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
                                         reference_gas_state=fill_state_gas_object)

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
            print ("Current guess for {0} = {1:.2f} m/s".format(shock_name, vs))

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

            v3g = unsteadily_expanding_state_gas_flow.finite_wave_dp(unsteadily_expanding_state.get_gas_state(), unsteadily_expanding_state.get_v(),
                                                                     'cplus', p3, unsteadily_expanded_gas_state, steps=steps)

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

            print("Current p{0} = {1:.2f} Pa, current p{2} = {3:.2f} Pa.".format(shocked_state_label_number, shocked_gas_state.p,
                                                                         unsteadily_expanded_state_label_number, unsteadily_expanded_gas_state.p))
            print("Current v{0}g = {1:.2f} m/s, current v{2}g = {3:.2f} m/s.".format(shocked_state_label_number, v2g,
                                                                             unsteadily_expanded_state_label_number, v3g))
            if abs((v2g - v3g) / v2g) > 0.001:
                print("Current (v{0}g - v{1}g) / v{0}g = {2:.6f}.".format(shocked_state_label_number, unsteadily_expanded_state_label_number, (v2g - v3g) / v2g))
            else:
                print("Current (v{0}g - v{1}g) / v{0}g = {2:.3e}.".format(shocked_state_label_number, unsteadily_expanded_state_label_number, (v2g - v3g) / v2g))

            return (v2g - v3g) / v2g

        print('-'*60)
        print("Calculating {0} shock speed ({1}).".format(self.tube_name, self.shock_name))
        print("{0} fill state is:".format(self.tube_name))
        print(self.fill_state)
        print("Unsteadily expanding entry state is:")
        print(self.unsteadily_expanding_state)

        # calculate the shock speed using a secant solver
        self.vs = secant(error_in_velocity_function, self.vs_guess_1, self.vs_guess_2, limits = self.vs_limits, tol=self.vs_tolerance)

        self.Ms = self.vs / self.fill_state.get_gas_state().a

        print ('-' * 60)
        print ("From secant solve: {0} = {1:.2f} m/s".format(self.shock_name, self.vs))

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
        print("Setting {0} shock speed ({1}) to {2:.2f} m/s".format(self.tube_name, self.shock_name, vs))

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
            print("Now that {0} is known, finding conditions at states {1} and {2}.".format(self.shock_name,
                                                                                            self.shocked_fill_state_name,
                                                                                            self.unsteadily_expanded_entrance_state_name))

            # now that we have the shock speed, we can get the final versions of the post-shock and unsteadily expanded states
            # (we don't use the values from the secant solver as we may want to change what we expand the unsteadily expanded state to)

            fill_state_gmodel = self.fill_state.get_gas_state().gmodel

            shocked_gas_state = GasState(fill_state_gmodel)
            fill_state_gas_flow = GasFlow(fill_state_gmodel)

            v2, v2g = fill_state_gas_flow.normal_shock(self.fill_state.get_gas_state(), self.vs, shocked_gas_state)

            # for the shocked state we can just use the fill state as the reference state...
            self.shocked_state = Facility_State(self.shocked_fill_state_name, shocked_gas_state, v2g,
                                                reference_gas_state=self.fill_state.get_gas_state())

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
                v3g = unsteadily_expanding_state_gas_flow.finite_wave_dp(self.unsteadily_expanding_state.get_gas_state(),
                                                                         self.unsteadily_expanding_state.get_v(),
                                                                         'cplus', p3, unsteadily_expanded_gas_state,
                                                                         steps=self.unsteady_expansion_steps)
            elif self.expand_to == 'shock_speed':
                # we use finite wave_dv and expand to the shock speed instead...
                v3g = unsteadily_expanding_state_gas_flow.finite_wave_dv(self.unsteadily_expanding_state.get_gas_state(),
                                                                         self.unsteadily_expanding_state.get_v(),
                                                                         'cplus', self.vs, unsteadily_expanded_gas_state,
                                                                         steps=self.unsteady_expansion_steps)

            # for this state, if the entrance state had a reference gas state, we can grab it as the reference state...

            if self.unsteadily_expanding_state.reference_gas_state:
                reference_gas_state = self.unsteadily_expanding_state.get_reference_gas_state()
            else:
                reference_gas_state = None
            self.unsteadily_expanded_state = Facility_State(self.unsteadily_expanded_entrance_state_name, unsteadily_expanded_gas_state, v3g,
                                                            reference_gas_state=reference_gas_state)

            for facility_state in [self.shocked_state, self.unsteadily_expanded_state]:
                print('-'*60)
                print(facility_state)
                if facility_state.get_gas_state().gmodel.type_str == 'CEAGas':
                    print('species in {0} at equilibrium (by mass):'.format(facility_state.get_state_name()))
                    print(facility_state.get_reduced_species_massf_dict())

                # add the stagnation enthalpy, if we can:
                if facility_state.reference_gas_state:
                    total_enthalpy = facility_state.get_total_enthalpy()

                    print ("The total enthalpy (Ht) at state {0} is {1:.2f} MJ/kg.".format(facility_state.get_state_name(),
                                                                                           total_enthalpy/1.0e6))
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
                if facility_state.get_gas_state().gmodel.type_str == 'CEAGas':
                    print('species in {0} at equilibrium (by mass):'.format(facility_state.get_state_name()))
                    print(facility_state.get_reduced_species_massf_dict())

                # add the stagnation enthalpy, if we can:
                if facility_state.reference_gas_state:
                    total_enthalpy = facility_state.get_total_enthalpy()

                    print ("The total enthalpy (Ht) at state {0} is {1:.2f} MJ/kg.".format(facility_state.get_state_name(),
                                                                                           total_enthalpy/1.0e6))


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
            return self.shocked_state
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

    def __init__(self, entrance_state_name, entrance_state, exit_state_name, area_ratio, nozzle_expansion_tolerance):
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

        # TO DO: maybe add some comments about what it is doing...

        print ('-'*60)
        print ("Starting steady expansion through the nozzle using an area ratio of {0}".format(self.area_ratio))

        # make the gas model object

        entrance_state_gmodel = self.entrance_state.get_gas_state().gmodel
        entrance_state_gas_flow_object = GasFlow(entrance_state_gmodel)

        exit_gas_state = GasState(entrance_state_gmodel)

        v_exit = entrance_state_gas_flow_object.steady_flow_with_area_change(self.entrance_state.get_gas_state(), self.entrance_state.get_v(),
                                                                             self.area_ratio, exit_gas_state,
                                                                             tol=nozzle_expansion_tolerance)

        # if the entrance state has a reference gas state, we can grab that as the exit state will have the same one.
        if self.entrance_state.reference_gas_state:
            reference_gas_state = self.entrance_state.get_reference_gas_state()
        else:
            reference_gas_state = None
        self.exit_state = Facility_State(self.exit_state_name, exit_gas_state, v_exit,
                                         reference_gas_state=reference_gas_state)

        print (self.exit_state)
        if self.exit_state.get_gas_state().gmodel.type_str == 'CEAGas':
            print('species in {0} at equilibrium (by mass):'.format(self.exit_state.get_state_name()))
            print(self.exit_state.get_reduced_species_massf_dict())

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

        # TO DO: work out how to get a perfect gas version, which was easy in PITOT

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
        self.post_normal_shock_state = Facility_State('{0}e'.format(self.test_section_post_shock_state_name),
                                                      post_normal_shock_gas_state, v10,
                                                      reference_gas_state=reference_gas_state)

        print (self.post_normal_shock_state)

        # TO DO: need to get the mole fractions working too that is normally what we want...
        if self.post_normal_shock_state.get_gas_state().gmodel.type_str == 'CEAGas':
            print('species in {0} at equilibrium (by mass):'.format(self.post_normal_shock_state.get_state_name()))
            print(self.post_normal_shock_state.get_reduced_species_massf_dict())

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

    def calculate_post_conical_shock_state(self, cone_half_angle_degrees):
        """
        Calculate the post-conical shock state for a given cone half angle (in degrees)

        :param cone_half_angle_degrees:
        :return:
        """

        self.cone_half_angle_degrees = cone_half_angle_degrees
        cone_half_angle = math.radians(self.cone_half_angle_degrees)

        print('-' * 60)
        print("Starting equilibrium conical shock calculation with a cone half angle of {0} degrees.".format(self.cone_half_angle_degrees))

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

        print("Shock angle over the cone is {0} degrees".format(self.conical_shock_half_angle_degrees))

        cone_half_angle_calculated, v10c = test_section_state_gas_flow_object.theta_cone(self.test_section_state.get_gas_state(),
                                                                                   self.test_section_state.get_v(),
                                                                                   beta,
                                                                                   post_conical_shock_gas_state)

        cone_half_angle_calculated_degrees = math.degrees(cone_half_angle_calculated)

        print ("The calculated cone half-angle should be the same as the specified one: {0} deg = {1} deg".format(cone_half_angle_calculated_degrees,
                                                                                                                  self.cone_half_angle_degrees))

        # if the entrance state has a reference gas state, we can grab that as the post-shock state will have the same one.
        if self.entrance_state.reference_gas_state:
            reference_gas_state = self.entrance_state.get_reference_gas_state()
        else:
            reference_gas_state = None

        self.post_conical_shock_state = Facility_State('{0}c'.format(self.test_section_post_shock_state_name),
                                                       post_conical_shock_gas_state, v10c,
                                                       reference_gas_state=reference_gas_state)

        print(self.post_conical_shock_state)

        if self.post_normal_shock_state.get_gas_state().gmodel.type_str == 'CEAGas':
            print('species in {0} at equilibrium (by mass):'.format(self.post_normal_shock_state.get_state_name()))
            print(self.post_normal_shock_state.get_reduced_species_massf_dict())

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
        print("Starting equilibrium wedge shock calculation with a wedge angle of {0} degrees.".format(self.wedge_angle_degrees))

        print("Test section freestream state is:")
        print(self.test_section_state)

        test_section_state_gmodel = self.test_section_state.get_gas_state().gmodel
        test_section_state_gas_flow_object = GasFlow(test_section_state_gmodel)

        post_wedge_shock_gas_state = GasState(test_section_state_gmodel)

        # start by getting the shock angle (beta)
        beta = test_section_state_gas_flow_object.beta_oblique(self.test_section_state.get_gas_state(),
                                                               self.test_section_state.get_v(),
                                                               wedge_angle)

        self.wedge_shock_angle_degrees = math.degrees(beta)

        print("Shock angle over the wedge is {0} degrees".format(self.wedge_shock_angle_degrees))

        wedge_angle_calculated, v10w = test_section_state_gas_flow_object.theta_oblique(self.test_section_state.get_gas_state(),
                                                                                        self.test_section_state.get_v(),
                                                                                        beta,
                                                                                        post_wedge_shock_gas_state)

        wedge_angle_calculated_degrees = math.degrees(wedge_angle_calculated)

        # TO DO: could add the check on this angle here which PITOT had...

        print ("The calculated wedge angle should be the same as the specified one: {0} deg = {1} deg".format(wedge_angle_calculated_degrees,
                                                                                                               self.wedge_angle_degrees))

        # if the entrance state has a reference gas state, we can grab that as the post-shock state will have the same one.
        if self.entrance_state.reference_gas_state:
            reference_gas_state = self.entrance_state.get_reference_gas_state()
        else:
            reference_gas_state = None

        self.post_wedge_shock_state = Facility_State('{0}w'.format(self.test_section_post_shock_state_name),
                                                     post_wedge_shock_gas_state, v10w,
                                                     reference_gas_state=reference_gas_state)

        print(self.post_wedge_shock_state)

        if self.post_wedge_shock_state.get_gas_state().gmodel.type_str == 'CEAGas':
            print('species in {0} at equilibrium (by mass):'.format(self.post_wedge_shock_state.get_state_name()))
            print(self.post_wedge_shock_state.get_reduced_species_massf_dict())

        return

    def get_facility_states(self):
        """
        This function returns a list of the unique facility states (i.e. it shouldn't return anything
        which would be returned using the function with the same name on another object) which this object contains,
        in an order which is relevant for the PITOT3 text output.

        This was mainly made for the output, but may have other uses.

        """

        facility_state_list = [self.post_normal_shock_state]

        if hasattr(self, 'post_conical_shock_state'):
            facility_state_list += [self.post_conical_shock_state]

        if hasattr(self, 'post_wedge_shock_state'):
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
    if pitot_p/1000.0 < 10000.0:
        output_line += "{0:<8.1f}".format(pitot_p/1000.0) # to get kPa
    else:
        output_line += "{0:<8.0f}".format(pitot_p/1000.0)  # to get kPa
    if p0/1.0e6 < 1000.0:
        output_line += "{0:<7.2f}".format(p0/1.0e6) # to get MPa
    elif 1000.0 <= p0/1.0e6 < 10000.0:
        output_line += "{0:<7.1f}".format(p0 / 1.0e6)  # to get MPa
    else:
        output_line += "{0:<7.0f}".format(p0 / 1.0e6)  # to get MPa
    if Ht == '-':
        output_line += "{0:<6}".format(Ht)
    else:
        if Ht/1.0e6 < 100.0:
            output_line += "{0:<7.2f}".format(Ht/1.0e6) # to get MJ/kg
        else:
            output_line += "{0:<7.1f}".format(Ht/1.0e6)  # to get MJ/kg
    if h == '-':
        output_line += "{0:<6}".format(h)
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

def pitot3_results_output(config_data, gas_path, object_dict):
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

    if 'output_filename' in config_data:
        if '.txt' not in config_data['output_filename']:
            output_file = open(config_data['output_filename'] + '.txt', "w")
        else:
            output_file = open(config_data['output_filename'], "w")
        output_list.append(output_file)

    for output_stream in output_list:
        # words starting with vowels causing issues below...
        if config_data['facility_type'] in ['expansion_tube']:
            print("PITOT3 Version {0} doing an {1} calculation.".format(config_data['VERSION_STRING'], config_data['facility_type']),
                  file=output_stream)
        else:
            print("PITOT3 Version {0} doing a {1} calculation.".format(config_data['VERSION_STRING'], config_data['facility_type']),
                  file=output_stream)
        print("Calculation mode is '{0}'.".format(config_data['mode']), file=output_stream)
        if config_data['facility']:
            print("Facility is '{0}'.".format(config_data['facility']), file=output_stream)

        driver = object_dict['driver']
        state4 = driver.get_driver_rupture_state()
        if config_data['driver_condition'] != 'custom':
            print("Driver condition is '{0}'. Driver gas model is {1}.".format(config_data['driver_condition'],
                                                                               state4.get_gas_state().gmodel.type_str),
                  file=output_stream)
        else:
            print("Using custom driver condition from the file {0}.".format(config_data['driver_condition_filename']),
                  file=output_stream)
            print("Driver gas model is {0}.".format(state4.get_gas_state().gmodel.type_str),
                  file=output_stream)

        if state4.get_gas_state().gmodel.type_str == 'CEAGas':
            print("Driver gas composition is {0} (by mass) ({1}).".format(state4.get_reduced_species_massf_dict(),
                                                                          state4.get_gamma_and_R_string()),
                  file=output_stream)
        else:
            print("Driver gas {0}.".format(state4.get_gamma_and_R_string()), file=output_stream)
        if 'nozzle' in object_dict:
            print("Nozzle area ratio is {0}.".format(object_dict['nozzle'].get_area_ratio()), file=output_stream)

        if 'secondary_driver' in object_dict:
            secondary_driver = object_dict['secondary_driver']

            secondary_driver_fill_state = secondary_driver.get_fill_state()

            if secondary_driver.get_fill_gas_model() != 'custom':
                print("Secondary driver gas ({0}) is {1}. Secondary driver gas model is {2}.".format(
                    secondary_driver_fill_state.get_state_name(),
                    secondary_driver.get_fill_gas_name(),
                    secondary_driver_fill_state.get_gas_state().gmodel.type_str),
                      file=output_stream)
            else:
                print('Using custom secondary driver gas from the file {0}.'.format(secondary_driver.get_fill_gas_filename()),
                      file=output_stream)
                print("Secondary driver gas model is {0}.".format(secondary_driver_fill_state.get_gas_state().gmodel.type_str),
                    file=output_stream)

            if secondary_driver_fill_state.get_gas_state().gmodel.type_str == 'CEAGas':
                print("Secondary driver gas composition is {0} by mass ({1}).".format(
                    secondary_driver_fill_state.get_reduced_species_massf_dict(),
                    secondary_driver_fill_state.get_gamma_and_R_string()),
                      file=output_stream)
            else:
                print("Secondary driver gas {0}.".format(secondary_driver_fill_state.get_gamma_and_R_string()),
                      file=output_stream)

        # we need the test gas here, which is the shock tube fill state...

        shock_tube = object_dict['shock_tube']

        shock_tube_fill_state = shock_tube.get_fill_state()

        if shock_tube.get_fill_gas_model() != 'custom':
            print("Test gas ({0}) is {1}. Test gas gas model is {2}.".format(shock_tube_fill_state.get_state_name(),
                                                                             shock_tube.get_fill_gas_name(),
                                                                             shock_tube_fill_state.get_gas_state().gmodel.type_str),
                  file=output_stream)
        else:
            print('Using custom test gas from the file {0}.'.format(shock_tube.get_fill_gas_filename()),
                  file=output_stream)
            print("Test gas gas model is {0}.".format(shock_tube_fill_state.get_gas_state().gmodel.type_str),
                  file=output_stream)

        if shock_tube_fill_state.get_gas_state().gmodel.type_str == 'CEAGas':

            print("Test gas composition is {0} by mass ({1})".format(
                shock_tube_fill_state.get_reduced_species_massf_dict(),
                shock_tube_fill_state.get_gamma_and_R_string()),
                  file=output_stream)
        else:
            print("Test gas {0}".format(shock_tube_fill_state.get_gamma_and_R_string()),
                  file=output_stream)

        # if we have an acceleration tube, we are an expansion tube...
        if 'acceleration_tube' in object_dict:
            acceleration_tube = object_dict['acceleration_tube']

            acceleration_tube_fill_state = acceleration_tube.get_fill_state()

            if acceleration_tube.get_fill_gas_model() != 'custom':
                print("Accelerator gas ({0}) is {1}. Accelerator gas gas model is {2}.".format(
                    acceleration_tube_fill_state.get_state_name(),
                    acceleration_tube.get_fill_gas_name(),
                    acceleration_tube_fill_state.get_gas_state().gmodel.type_str),
                      file=output_stream)
            else:
                print('Using custom accelerator gas from the file {0}.'.format(acceleration_tube.get_fill_gas_filename()),
                      file=output_stream)
                print("Accelerator gas gas model is {0}.".format(acceleration_tube_fill_state.get_gas_state().gmodel.type_str),
                    file=output_stream)

            if acceleration_tube_fill_state.get_gas_state().gmodel.type_str == 'CEAGas':
                print("Accelerator gas composition is {0} by mass ({1}).".format(
                    acceleration_tube_fill_state.get_reduced_species_massf_dict(),
                    acceleration_tube_fill_state.get_gamma_and_R_string()),
                      file=output_stream)
            else:
                print("Accelerator gas {0}.".format(acceleration_tube_fill_state.get_gamma_and_R_string()),
                      file=output_stream)

        if 'secondary_driver' in locals():
            print('vsd = {0:.2f} m/s, Msd = {1:.2f}'.format(secondary_driver.get_shock_speed(),
                                                            secondary_driver.get_shock_Mach_number()),
                  file=output_stream)

        shock_speed_output = 'vs1 = {0:.2f} m/s, Ms1 = {1:.2f}'.format(shock_tube.get_shock_speed(),
                                                                       shock_tube.get_shock_Mach_number())

        if 'acceleration_tube' in locals():
            shock_speed_output += ', vs2 = {0:.2f} m/s, Ms2 = {1:.2f}'.format(acceleration_tube.get_shock_speed(),
                                                                              acceleration_tube.get_shock_Mach_number())

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
                tube_name_reduced = '{0}{1}'.format(tube_name.split('_')[0][0], tube_name.split('_')[1][0])
                vr, Mr = object_dict[diaphragm].get_vr_Mr()

                print("NOTE: a user specified reflected shock was done at the end of the {0}.".format(tube_name),
                      file=output_stream)

                print("vr-{0} = {1:.2f} m/s, Mr-{0} = {2:.2f}".format(tube_name_reduced, vr, Mr),
                      file=output_stream)

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

        # now do some final extra outputs at the bottom...
        # start by pulling out our freestream and post-shock test section states which we will need later...
        test_section = object_dict['test_section']

        freestream_state = test_section.get_entrance_state()
        test_section_post_normal_shock_state = test_section.get_post_normal_shock_state()

        freestream_total_temperature = freestream_state.get_total_temperature()

        freestream_flight_equivalent_velocity = freestream_state.get_flight_equivalent_velocity()

        print("The freestream ({0}) total temperature (Tt) is {1:.2f} K.".format(freestream_state.get_state_name(),
                                                                                 freestream_total_temperature),
              file=output_stream)
        print("The freestream ({0}) flight equivalent velocity (Ue) is {1:.2f} m/s.".format(
            freestream_state.get_state_name(),
            freestream_flight_equivalent_velocity),
              file=output_stream)

        if acceleration_tube in locals():
            if acceleration_tube.tube_length:
                basic_test_time = expansion_tube_test_time_calculator(acceleration_tube)
                print("Basic test time = {0:.2f} microseconds.".format(basic_test_time * 1.0e6), file=output_stream)

        if hasattr(test_section, 'post_wedge_shock_state'):
            print("Wedge angle was {0:.2f} degrees. Wedge shock angle (beta) was found to be {1:.2f} degrees.".format(
                test_section.wedge_angle_degrees,
                test_section.wedge_shock_angle_degrees),
                  file=output_stream)

        print("Freestream ({0}) unit Reynolds number is {1:.2f} /m (related mu is {2:.2e} Pa.s)." \
              .format(freestream_state.get_state_name(), freestream_state.get_unit_Reynolds_number(),
                      freestream_state.get_mu()),
              file=output_stream)

        print("Post normal shock equilibrium ({0}) unit Reynolds number is {1:.2f} /m (related mu is {2:.2e} Pa.s)." \
              .format(test_section_post_normal_shock_state.get_state_name(),
                      test_section_post_normal_shock_state.get_unit_Reynolds_number(),
                      test_section_post_normal_shock_state.get_mu()),
              file=output_stream)

        if freestream_state.get_gas_state().gmodel.type_str == 'CEAGas':
            print("Species in the freestream state ({0}) at equilibrium (by mass):".format(
                freestream_state.get_state_name()),
                  file=output_stream)
            print(freestream_state.get_reduced_species_massf_dict(), file=output_stream)

        if test_section.get_post_normal_shock_state().get_gas_state().gmodel.type_str == 'CEAGas':
            print("Species in the shock layer at equilibrium ({0}) (by mass):".format(
                test_section.get_post_normal_shock_state().get_state_name()),
                  file=output_stream)
            print(test_section.get_post_normal_shock_state().get_reduced_species_massf_dict(), file=output_stream)

    return
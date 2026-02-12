# reflected_shock_tunnel.py
"""
Author
------
PA Jacobs
Institute of Aerodynamics and Flow Technology
The German Aerospace Center, Goettingen.
                 &
School of Mechancial & Mining Engineering,
The University of Queensland

Versions
--------
24-Dec-2002 PJ: First code (estcj.py).
2010 PJ : ported to run with Rowan's cea2_gas module.
2011 PJ : Added isentropic expansions so that we now have
    a full replacement for stn.f
01-June-2011 Luke Doherty: Separated the code which writes an output
    file into its own function to allow for better integration with nenzfr.py
30-June-2011 Luke Doherty: Decreased the starting guess for secant
    when solving for the exit flow
22-July-2011 Luke Doherty: Added stnp option which allows us to expand
    to a nominated pitot-to-supply pressure ratio. The calculated pitot
    pressure and pitot-to-supply pressure ratio are included in the values
    printed out for the nozzle exit
24-Feb-2012 PJ: update to use the new cea2_gas.py arrangement.
31-Dec-2013 PJ: added libgas_gas.py option.
14-Jan-2014 PJ: included ideal gas option.
06-Jun-2020 PJ: Updated to use the loadable library from Eilmer4
"""


from gdtk.gas import GasModel, GasState, GasFlow
from gdtk.numeric.zero_solvers import secant


def calculate_states(gasModel, flow, p1, T1, massf, Vs, pe, pp_on_pe, area_ratio,
                     task, print_status=True):
    """
    Runs the reflected-shock-tube calculation from initial fill conditions
    observed shock speed and equilibrium pressure.

    Input
    -----
    gasModel: gas model object
    flow: flow functions object
    p1: fill pressure of gas initially filling shock tube
    T1: fill temperature of gas initially filling shock tube
    massf: list of mass fractions
    Vs: observed incident shock speed
    pe: observed pressure once shock-reflected region reaches equilibrium
    pp_on_pe: specify this ratio if we want the supersonic nozzle expansion to
        terminate at a particular Pitot pressure
    area_ratio: specify this ratio if we want the supersonic nozzle expansion
        to proceed to a particular quasi-one-dimensional area ratio.
    task: one of 'ishock', 'st', 'stn', 'stnp'
    print_status: if True, the start of each stage of the computation is noted.

    Returns a dictionary of results, including gas states and velocities.
    """
    MY_DEBUG = False  # some detailed data is output to help debugging
    #
    if print_status: print('Write pre-shock condition.')
    state1 = GasState(gasModel)
    state1.p = p1; state1.T = T1; state1.massf = massf
    state1.update_thermo_from_pT()
    H1 = state1.internal_energy + state1.p/state1.rho
    result = {'state1':state1, 'H1':H1}
    #
    if print_status: print('Start incident-shock calculation.')
    state2 = GasState(gasModel)
    (V2,Vg) = flow.normal_shock_1(state1, Vs, state2)
    result['state2'] = state2
    result['V2'] = V2
    result['Vg'] = Vg
    #
    if task == 'ishock':
        # We want post-incident-shock conditions only.
        return result
    #
    if print_status: print('Start reflected-shock calculation.')
    state5 = GasState(gasModel)
    Vr = flow.reflected_shock(state2, Vg, state5)
    result['state5'] = state5
    result['Vr'] = Vr
    #
    if print_status: print('Start calculation of isentropic relaxation.')
    state5s = GasState(gasModel); state5s.copy_values(state5)
    # entropy is set, then pressure is relaxed via an isentropic process
    state5s.p = state5.p if pe == None else pe
    if MY_DEBUG: print("state5.entropy=", state5.entropy)
    state5s.update_thermo_from_ps(state5.entropy)
    result['state5s'] = state5s
    H5s = state5s.internal_energy + state5s.p/state5s.rho # stagnation enthalpy
    result['H5s'] = H5s
    #
    if task in ['stn','stnp']:
        if print_status: print('Start isentropic relaxation to throat (Mach 1)')
        def error_at_throat(x, s5s=state5s):
            "Returns Mach number error as pressure is changed."
            state = GasState(gasModel)
            V = flow.expand_from_stagnation(s5s, x, state)
            return (V/state.a) - 1.0
        try:
            x6 = secant(error_at_throat, 0.95, 0.90, tol=1.0e-4)
        except RuntimeError:
            print("Failed to find throat conditions iteratively.")
            x6 = 1.0
        state6 = GasState(gasModel)
        V6 = flow.expand_from_stagnation(state5s, x6, state6)
        mflux6 = state6.rho * V6  # mass flux per unit area, at throat
        result['state6'] = state6
        result['V6'] = V6
        result['mflux6'] = mflux6
        #
        if task == 'stn':
            if print_status: print('Start isentropic relaxation to nozzle exit of given area.')
            # The mass flux going through the nozzle exit has to be the same
            # as that going through the nozzle throat.
            def error_at_exit(x, s5s=state5s, s6=state6, mflux_throat=mflux6,
                              area_ratio=area_ratio):
                "Returns mass_flux error as pressure is changed."
                state = GasState(gasModel)
                V = flow.expand_from_stagnation(s5s, x, state)
                mflux = state.rho * V * area_ratio
                if MY_DEBUG:
                    print("x=", x, "p=", state.p, "T=", state.T, "V=", V, \
                          "mflux=", mflux, "mflux_throat=", mflux_throat)
                return (mflux-mflux_throat)/mflux_throat
            # It appears that we need a pretty good starting guess for the pressure ratio.
            # Maybe a low value is OK.
            try:
                x7 = secant(error_at_exit, 0.001*x6, 0.00005*x6, tol=1.0e-4,
                            limits=[1.0/state5s.p,1.0])
            except RuntimeError:
                print("Failed to find exit conditions iteratively.")
                x7 = x6
            state7 = GasState(gasModel)
            V7 = flow.expand_from_stagnation(state5s, x7, state7)
            mflux7 = state7.rho * V7 * area_ratio
            result['area_ratio'] = area_ratio
            state7_pitot = GasState(gasModel)
            flow.pitot_condition(state7, V7, state7_pitot)
            result['state7'] = state7
            result['V7'] = V7
            result['mflux7'] = mflux7
            result['pitot7'] = state7_pitot.p
        elif task == 'stnp':
            if print_status: print('Start isentropic relaxation to nozzle exit pitot pressure.')
            # The exit pitot pressure has to be the same as that measured
            def error_at_exit(x, s5s=state5s, s6=state6, pp_pe=pp_on_pe):
                "Returns pitot pressure error as static pressure is changed."
                state1 = GasState(gasModel)
                V = flow.expand_from_stagnation(s5s, x, state1)
                state2 = GasState(gasModel)
                flow.pitot_condition(state1, V, state2)
                if MY_DEBUG:
                    print("x=", x, "pitot_to_supply=", state2.p/s5s.p, \
                          "relative error=", (state2.p/s5s.p - pp_pe)/pp_pe)
                return (state2.p/s5s.p - pp_pe)/pp_pe
            # We need a low starting guess for the pressure ratio.
            #x7 = secant(error_at_exit, 0.001*x6, 0.00005*x6, tol=1.0e-4)
            # Changed the tolerance on 25/07/2011 in order to get the M8 nozzle to work (shot 10803)
            try:
                x7 = secant(error_at_exit, 0.001*x6, 0.00005*x6, tol=2.0e-4,
                            limits=[1.0/state5s.p,1.0])
            except RuntimeError:
                print("Failed to find exit conditions iteratively.")
                x7 = x6
            state7 = GasState(gasModel)
            V7 = flow.expand_from_stagnation(state5s, x7, state7)
            result['area_ratio'] = mflux6/(state7.rho * V7)
            state7_pitot = GasState(gasModel)
            flow.pitot_condition(state7, V7, state7_pitot)
            result['state7'] = state7
            result['V7'] = V7
            result['mflux7'] = mflux6
            result['pitot7'] = state7_pitot.p
            if MY_DEBUG: print("area_ratio=", area_ratio, "pitot7=", state7_pitot.p)
    #
    if print_status: print('Done with reflected shock tube calculation.')
    return result

#! /usr/bin/env python3
"""
estcn.py: Equilibrium Shock Tube Conditions, with Nozzle

This program can be used to estimate flow conditions
for shock-processed flows typical of high-performance
shock-tunnels and expansion tubes.

The gas is assumed to remain in thermochemical equilibrium
and the flow processing is done in decoupled quasi-one-dimensional
wave processes such as shock waves and expansion fans.
For the reflected shock tunnel, this means that the initial,
quiescent test gas is first processed by the incident shock
and subsequently by the reflected shock.
The incident shock sets the inflow conditions for the reflected shock
but there is no further interaction.

The program can do a number of calculations:
* flow in a reflected shock tube with or without a nozzle
* pitot pressure from free-stream flow condition
* stagnation (total) condition from free-stream condition
* code surface condition from free-stream condition

When run as an application, this program takes its input as
command line arguments, performs the requested calculations and outputs
the gas-state results.
To see what specific inputs are required, start the program as::

$ estcn --help

Which particular input parameters you need to supply depends on the
chosen task, however, a typical flow condition for the T4 shock tunnel
with the Mach 4 nozzle may be computed using::

$ estcn --task=stn --gas=cea-lut-air.lua \
        --T1=300 --p1=125.0e3 --Vs=2414 --pe=34.37e6 --ar=27.0

The full output is a bit too much to include here, but you should see that
this condition has an enthalpy of 5.43 MJ/kg and the nozzle-exit condition
has a pressure of 93.6 kPa and a static temperature of 1284 degrees K,
with a flow speed of 2.95 km/s.


Some History
------------
Since 1968, we have been using the ESTC code by Malcolm McIntosh
to compute the conditions in the end of the reflected shock tubes
T1--T5 and HEG.  There are a number of problems in using the ESTC
code, including uncertainty in updating the chemistry coefficients.
This program, ESTCn, moves away from the old chemistry model
by making use of the gas models built into the Eilmer4 flow code.

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
06-Jun-2020 PJ: Renamed to estcn, updated to use the loadable library from Eilmer4
11-Jun-2020 PJ: State calculations for the reflected-shock-tunnel moved to a separate module.
    All that is left here is the code for gathering input and presenting results.
"""

VERSION_STRING = "11-Jun-2020"

import sys, os, math
from eilmer.gas import GasModel, GasState, GasFlow
from eilmer.reflected_shock_tunnel import calculate_states


def str_with_units(gs, lead_str='  '):
    """
    Returns a string describing the gas state.  Quantities and units.
    """
    result = lead_str+'p: %g Pa, T: %g K, rho: %g kg/m**3, u: %g J/kg, h: %g J/kg' % \
             (gs.p, gs.T, gs.rho, gs.internal_energy, gs.enthalpy)
    result += '\n'+lead_str+'R: %g J/(kg.K), gam: %g, Cp: %g J/(kg.K), a: %g m/s, s: %g J/(kg.K)' % \
              (gs.R, gs.gamma, gs.Cp, gs.a, gs.entropy)
    if gs.gmodel.type_str == "CEAGas":
        result += '\n'+lead_str+('CEA-massf: %s' % gs.ceaSavedData['massf'])
    return result


def main():
    """
    The application gets information from the command options,
    does some calculation (depending on the specified task)
    and writes the results to the console or a file.
    """
    import optparse
    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('--task', dest='task', default='st',
                  choices=['st', 'stn', 'stnp', 'ishock', 'total', 'pitot', 'cone'],
                  help=("particular calculation to make: "
                        "st = reflected shock tube; "
                        "stn = reflected shock tube with nozzle; "
                        "stnp = reflected shock tube with nozzle expanded to pitot; "
                        "ishock = incident shock only; "
                        "total = free-stream to total condition; "
                        "pitot = free-stream to Pitot condition; "
                        "cone = free-stream to Taylor-Maccoll cone flow"))
    op.add_option('--gas-model', dest='gasFileName', default='cea-lut-air.lua',
                  help=("name of gas model file"))
    op.add_option('--p1', dest='p1', type='float', default=None,
                  help=("shock tube fill pressure or static pressure, in Pa"))
    op.add_option('--T1', dest='T1', type='float', default=None,
                  help=("shock tube fill temperature, in degrees K"))
    op.add_option('--molef', dest='molef_str', type='string', default=None,
                  help=("mole fractions of gas model species,"
                        " provided as a comma-separated list of float values"))
    op.add_option('--massf', dest='massf_str', type='string', default=None,
                  help=("mass fractions of gas model species,"
                        " provided as a comma-separated list of float values"))
    op.add_option('--V1', dest='V1', type='float', default=None,
                  help=("initial speed of gas in lab frame [default: %default], in m/s"))
    op.add_option('--Vs', dest='Vs', type='float', default=None,
                  help=("incident shock speed, in m/s"))
    op.add_option('--pe', dest='pe', type='float', default=None,
                  help=("equilibrium pressure (after shock reflection), in Pa"))
    op.add_option('--pp_on_pe', dest='pp_on_pe', type='float', default=None,
                  help=("nozzle supply to exit pitot pressure ratio"))
    op.add_option('--ar', dest='area_ratio', type='float', default=None,
                  help=("exit-to-throat area ratio of the nozzle"))
    op.add_option('--sigma-deg', dest='cone_half_angle_deg', type='float', default=None,
                  help=("half-angle of the cone, in degrees"))
    op.add_option('--ofn', dest='outFileName', default=None,
                  help="name of file in which to accumulate output."
                      " file name will be: outFileName-estcn.dat"
                      " (Note that output defaults to stdout.)")
    opt, args = op.parse_args()
    #
    task = opt.task
    gasFileName = opt.gasFileName
    p1 = opt.p1
    T1 = opt.T1
    molef_str = opt.molef_str
    massf_str = opt.massf_str
    V1 = opt.V1
    Vs = opt.Vs
    pe = opt.pe
    pp_on_pe = opt.pp_on_pe
    area_ratio = opt.area_ratio
    cone_half_angle_deg = opt.cone_half_angle_deg
    outFileName = opt.outFileName
    #
    bad_input = False
    if not os.path.exists(gasFileName):
        print("Gas model file "+gasFileName+" seems not to exist.")
        bad_input = True
    if p1 is None:
        print("Need to supply a float value for p1.")
        bad_input = True
    if T1 is None:
        print("Need to supply a float value for T1.")
        bad_input = True
    if Vs is None and task in ['stn', 'stnp', 'st', 'ishock']:
        print("Need to supply a float value for Vs.")
        bad_input = True
    if V1 is None and task in ['pitot', 'total', 'cone']:
        print("Need to supply a free-stream velocity.")
        bad_input = True
    if cone_half_angle_deg is None and task in ['cone',]:
        print("Need to supply a cone half-angle (in degrees).")
        bad_input = True
    if pe is None and task in ['stn', 'stnp', 'st']:
        print("Need to supply a float value for pe.")
        bad_input = True
    if pp_on_pe is None and task in ['stnp']:
        print("Need to supply a float value for pp_on_pe.")
        bad_input = True
    if area_ratio is None and task in ['stn']:
        print("Need to supply a float value for ar=area_ratio.")
        bad_input = True
    if bad_input:
        return -2
    #
    fout = open(outFileName+'-estcn.dat','w') if outFileName else sys.stdout
    gasModel = GasModel(gasFileName)
    flow = GasFlow(gasModel)
    #
    massf = [1.0,] # default for nspecies == 1
    if gasModel.n_species > 1:
        if molef_str:
            molef = [float(s) for s in molef_str.split(',')]
            if abs(sum(molef)-1.0) > 1.0e-6:
                print("Mole fractions do not sum to 1.0; molef=", molef)
                return -2
            massf = gasModel.molef2massf(molef)
        elif massf_str:
            massf = [float(s) for s in massf_str.split(',')]
            if abs(sum(massf)-1.0) > 1.0e-6:
                print("Mass fractions do not sum to 1.0; massf=", massf)
                return -2
        else:
            print("Gas model requires species fractions but none supplied.")
            return -2
    #
    if task in ['st', 'stn', 'stnp', 'ishock']:
        # Some form of shock processing
        fout.write('Input parameters:\n')
        fout.write('  gasFileName is %s, p1: %g Pa, T1: %g K, massf=%s, Vs: %g m/s\n'
                   % (gasFileName, p1, T1, massf, Vs) )
        result = calculate_states(gasModel, flow, p1, T1, massf, Vs, pe, pp_on_pe, area_ratio, task=task)
        fout.write('State 1: pre-shock condition\n')
        fout.write(str_with_units(result['state1'])+'\n')
        fout.write('State 2: post-shock condition.\n')
        fout.write(str_with_units(result['state2'])+'\n')
        fout.write('  V2: %g m/s, Vg: %g m/s\n' % (result['V2'],result['Vg']) )
        if task in ['st', 'stn', 'stnp']:
            # Reflected-shock and, maybe, more
            fout.write('State 5: reflected-shock condition.\n')
            fout.write(str_with_units(result['state5'])+'\n')
            fout.write('  Vr: %g m/s\n' % (result['Vr'],) )
            fout.write('State 5s: equilibrium condition (relaxation to pe)\n')
            fout.write(str_with_units(result['state5s'])+'\n')
            fout.write('Enthalpy difference (H5s - H1): %g J/kg\n' %
                       ((result['H5s'] - result['H1']),) )
            if task in ['stn','stnp']:
                # Shock tube plus nozzle, expand gas isentropically to nozzle exit
                fout.write('State 6: Nozzle-throat condition (relaxation to M=1)\n')
                fout.write(str_with_units(result['state6'])+'\n')
                fout.write('  V6: %g m/s, M6: %g, mflux6: %g kg/s/m**2\n' %
                           (result['V6'], result['V6']/result['state6'].a, result['mflux6'],) )
                fout.write('State 7: Nozzle-exit condition (relaxation to correct mass flux)\n')
                fout.write(str_with_units(result['state7'])+'\n')
                fout.write('  V7: %g m/s, M7: %g, mflux7: %g kg/s/m**2, area_ratio: %g, pitot: %g Pa\n' %
                           (result['V7'], result['V7']/result['state7'].a, result['mflux7'],
                            result['area_ratio'], result['pitot7'],) )
                fout.write('  pitot7_on_p5s: %g\n' % (result['pitot7']/result['state5s'].p,) )
    elif task in ['total', 'TOTAL', 'Total']:
        # Isentropic total-pressure condition from free stream
        fout.write('Input parameters:\n')
        fout.write('  gasFileName is %s, p1: %g Pa, T1: %g K, massf: %s, V1: %g m/s\n'
                   % (gasFileName, p1, T1, massf, V1) )
        state1 = GasState(gasModel)
        state1.p = p1; state1.T = T1; state1.massf = massf
        state1.update_thermo_from_pT()
        state0 = GasState(gasModel); state0.copy_values(state1)
        flow.total_condition(state1, V1, state0)
        fout.write('Total condition:\n')
        fout.write(str_with_units(state0)+'\n')
    elif task in ['pitot', 'PITOT', 'Pitot']:
        # Pitot condition from free stream
        fout.write('Input parameters:\n')
        fout.write('  gasFileName is %s, p1: %g Pa, T1: %g K, massf: %s V1: %g m/s\n'
                   % (gasFileName, p1, T1, massf, V1) )
        state1 = GasState(gasModel)
        state1.p = p1; state1.T = T1; state1.massf = massf
        state1.update_thermo_from_pT()
        state0 = GasState(gasModel); state0.copy_values(state1)
        flow.pitot_condition(state1, V1, state0)
        fout.write('Pitot condition:\n')
        fout.write(str_with_units(state0)+'\n')
    elif task in ['cone', 'CONE', 'Cone']:
        # Conical shock processing from free stream
        fout.write('Input parameters:\n')
        fout.write('  gasFileName is %s, p1: %g Pa, T1: %g K, massf: %s, V1: %g m/s, sigma: %g degrees\n'
                   % (gasFileName, p1, T1, massf, V1, cone_half_angle_deg) )
        state1 = GasState(gasModel)
        state1.p = p1; state1.T = T1; state1.massf = massf
        state1.update_thermo_from_pT()
        fout.write('Free-stream condition:\n')
        fout.write(str_with_units(state1)+'\n')
        cone_half_angle_rad = math.radians(cone_half_angle_deg)
        beta_rad = flow.beta_cone(state1, V1, cone_half_angle_rad)
        state2 = GasState(gasModel); state2.copy_values(state1)
        theta_c, V_cone_surface = flow.theta_cone(state1, V1, beta_rad, state2)
        assert abs(theta_c - cone_half_angle_rad) < 0.001
        fout.write('Shock angle: %g (rad), %g (deg)\n' % (beta_rad, math.degrees(beta_rad)))
        fout.write('Cone-surface velocity: %g m/s\n' % (V_cone_surface,))
        fout.write('Cone-surface condition:\n')
        fout.write(str_with_units(state2)+'\n')
    #
    if outFileName: fout.close()
    return 0


if __name__ == '__main__':
    print("Equilibrium Shock Tube Conditions, with Nozzle")
    print("  Version:", VERSION_STRING)
    if len(sys.argv) <= 1:
        print("  To see some useful hints, invoke this program with option --help.")
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)

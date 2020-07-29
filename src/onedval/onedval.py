#!/usr/bin/env python
"""
Python program to compute one-dimensionalized quantities from 2D data.

This program is supplied with two inputs by the user:
1) a previously prepared Tecplot slice of data in ASCII block format; and
2) a config file to coordinate the calculation.
The program will then compute one-dimensionalized quantities from the
the 2D surface based on the methods selected by the user.

.. Author: Rowan J. Gollan

.. Versions:
   12-Feb-2013  Initial coding.
"""

import sys
from math import *
from gaspy import *
from cell import create_cells_from_slice, create_cells_from_line, area
from prop_avg import *
from copy import copy

output_formats = ['verbose', 'as_data_file']

default_var_map = {'x':'x', 'y':'y', 'z':'z', 'u':'u', 'v':'v', 'w':'w',
                   'rho':'rho', 'p':'p', 'T':'T', 'M':'M', 'h0':'h0'}

def print_usage():
    print ""
    print "Usage: onedval CONFIG INPUT_FILE(S)"
    print ""
    print "where:"
    print "CONFIG -- name of config file to control calculation"
    print "INPUT_FILE(S)  -- a list of one or more Tecplot files with slice data"
    print ""

# Global variable to hold gas model pointer
gmodel = None

pretty_var_names = {'rho':'density (kg/m^3)',
                    'p':'pressure (Pa)',
                    'p0':'total pressure (Pa)',
                    'T':'temperature (K)',
                    'T0':'total temperature (K)',
                    'M':'Mach number',
                    'u':'U velocity (m/s)',
                    'v':'V velocity (m/s)',
                    'w':'W velocity (m/s)'}

def pretty_print_props(f, props, species, outputs):
    for o in outputs:
        if o in pretty_var_names:
            f.write("%s\n" % pretty_var_names[o])
            f.write("%s = %.6e\n" % (o, props[o]))
        elif o in species:
            f.write("mass fraction of %s : %.6e\n" % (o, props[o]))
        else:
            f.write("%s : %.6e\n" % (o, props[o]))
        f.write("\n")
    return

short_one_d_av_names = {'area-weighted':'a-w',
                        'mass-flux-weighted':'m-w',
                        'flux-conserved':'f-c'}
units = {'rho':'kg/m^3',
         'p':'Pa',
         'p0':'Pa',
         'p0_rel':'Pa',
         'T':'K',
         'T0':'K',
         'T0_rel':'K',
         'a':'m/s',
         'R':'J/(kg.K)',
         'Cp':'J/(kg.K)',
         'Cv':'J/(kg.K)',
         'gamma':'-',
         'h':'J/kg',
         'h0':'J/kg',
         'h0_rel':'J/kg',
         'M':'-',
         'M_rel':'-',
         's':'J/(kg.K)',
         's_rel':'J/(kg.K)',
         'u':'m/s',
         'v':'m/s',
         'w':'m/s',
         'mass flux':'kg/s',
         'momentum flux':'kg.m/s^2',
         'energy flux':'W',
         'TKE':'W',
         'hr': 'J/kg'}

column_names = {'mass flux' : 'mass-flux',
                'momentum flux' : 'mom.-flux',
                'energy flux' : 'energy-flux' }

def data_file_header(f, one_d_av, one_d_out, int_out, species):
    f.write("# x[m]     y[m]     z[m]    area[m^2]     ")
    
    for av in one_d_av:
        for o in one_d_out:
            f.write("%s[%s]:%s    " % (o, units[o], short_one_d_av_names[av]))
    
    for i_out in int_out:
        if isinstance(i_out, tuple):
            f.write("%s    " % i_out[0])
        elif i_out == 'species mass flux':
            for sp in species:
                f.write("%s-flux[kg/s]    " % sp)
        else:
            f.write("%s[%s]    " % (column_names[i_out], units[i_out]))
    f.write("\n")
    return

def data_file_row(f, pos, area, phis, int_quants, one_d_av, one_d_out, int_out, species):
    f.write("%20.12e %20.12e %20.12e %20.12e " % (pos.x, pos.y, pos.z, area))
    for av in one_d_av:
        for o in one_d_out:
            f.write("%20.12e " % (phis[av][o]))
    
    for i_out in int_out:
        if isinstance(i_out, tuple):
            f.write("%20.12e " % int_quants[i_out[0]])
        elif i_out == 'species mass flux':
            for sp in species:
                f.write("%20.12e " % int_quants['species mass flux'][sp])
        else:
            f.write("%20.12e  " % (int_quants[i_out]))
    f.write("\n")


def main():
    """
    Top-level function for the onedval program.
    """
    print "onedval: A program to compute integrated and one-dimensionalised quantities."
    print "onedval: Beginning."
    # 0. Gather command-line info
    if len(sys.argv) < 3:
        print "At least two arguments are required."
        print_usage()
        sys.exit(1)

    config_file = sys.argv[1]
    
    # 1. Gather info from config file
    # Set some defaults.
    # If set to 'None', we expect to find something from the user.
    cfg = {}
    try:
        execfile(config_file, globals(), cfg)
    except IOError:
        print "There was a problem reading the config file: ", config_file
        print "Check that is conforms to Python syntax."
        print "Bailing out!"
        sys.exit(1)
    
    print "onedval: Setting up gas model"
    # 1a. Setup gas model
    if 'species' in cfg and 'gmodel_file' in cfg:
        print "There is a problem in the config file: ", config_file
        print "Both 'species' and 'gmodel_file' have been specified."
        print "You can only specify one of these keywords."
        print "Bailing out!"
        sys.exit(1)
    # Look for species.
    gmodel = None
    if 'species' in cfg:
        create_gas_file("thermally perfect gas", cfg['species'], "gas-model.lua")
        gmodel = create_gas_model("gas-model.lua")
        nsp = gmodel.get_number_of_species()
    # Else, test for gmodel_file
    if 'gmodel_file' in cfg:
        gmodel = create_gas_model(cfg['gmodel_file'])
        nsp = gmodel.get_number_of_species()
        cfg['species'] = []
        for isp in range(nsp):
            cfg['species'].append(gmodel.species_name(isp))

    if gmodel is None:
        print "There is a problem in the config file: ", config_file
        print "It appears that neither the 'species' nor the 'gmodel_file' have been specified."
        print "You must specify one of these keywords."
        print "Bailing out!"
        sys.exit(1)

    print "onedval: Checking over user inputs"
    # 1b. Check for variable map
    if not 'variable_map' in cfg:
        print "No 'variable_map' was set so the following default is used:"
        print default_var_map
        cfg['variable_map'] = default_var_map
    
    # 1c. Look for one_d_averages methods
    if not 'one_d_averages' in cfg:
        print "No list of 'one_d_averages' was set."
        print "The default method of 'flux-conserved' will be used."
        cfg['one_d_averages'] = ['flux-conserved']

    # 1d. Look for grid_scale
    if not 'grid_scale' in cfg:
        print "No 'grid_scale' was set."
        print "The default value of 1.0 will be used."
        cfg['grid_scale'] = 1.0

    # 1e. Look for one_d_outpus
    if not 'one_d_outputs' in cfg:
        cfg['one_d_outputs'] = []

    # 1f. Look for integrated_outputs
    if not 'integrated_outputs' in cfg:
        cfg['integrated_outputs'] = []

    # 1g. Looking for output options
    if not 'output_file' in cfg:
        print "No 'output_file' was set."
        print "An output file name must be set by the user."
        print "Bailing out!"
        sys.exit(1)
    
    if not 'output_format' in cfg:
        print "No 'output_format' was set."
        print "An output format must be set by the user."
        print "Bailing out!"
        sys.exit(1)

    if not cfg['output_format'] in output_formats:
        print "The selected output format: ", cfg['output_format']
        print "is not one of the available options. The available options are:"
        for o in output_formats:
            print "'%s'" % o
        print "Bailing out!"
        sys.exit(1)
        
    # 1h. Look for what to do with "bad" slices
    if not 'skip_bad_slices' in cfg:
        print "No 'skip_bad_slices' was set."
        print "The default value of True will be used."
        print "If bad data slices are encountered, the user"
        print "will be warned but the program will continue"
        print "on to the next slice and ignore the 'bad' slice."
        cfg['skip_bad_slices'] = True

    # 1i. Look for what kind of geometry
    if not 'geometry_type' in cfg:
        print "No 'geometry_type' was set."
        print "The default type of 3D will be used."
        cfg['geometry_type'] = '3D'

    # 2. Read data from slices and process
    print "onedval: Reading in data from slice(s)"
    f = open(cfg['output_file'], 'w')
    phis = {}
    int_quants = {}
    if cfg['output_format'] == 'as_data_file':
        data_file_header(f, cfg['one_d_averages'], cfg['one_d_outputs'], cfg['integrated_outputs'], cfg['species'])

    
    for slice_file in sys.argv[2:]:
        result = 'success'
        print "onedval: Creating cells from slice: ", slice_file
        if cfg['geometry_type'] == 'axi':
            cells = create_cells_from_line(slice_file, cfg['variable_map'], cfg['grid_scale'])
        else:
            cells = create_cells_from_slice(slice_file, cfg['variable_map'], cfg['grid_scale'])
        print "Total number of cells created from slice: ", len(cells)
        print "Make all cell normals consistent with very first cell."
        for c in cells[1:]:
             c._normal = cells[0]._normal
        # 2a. apply filtering if required
        if 'filter_function' in cfg:
            print "Using filter function to remove unwanted or unimportant cells."
            cells = filter(cfg['filter_function'], cells)
            print "Number of cells after filtering: ", len(cells)

        # 3. Do some work.
        print "onedval: Doing the requested calculations"
        if cfg['output_format'] == 'verbose':
            f.write("------------------- onedval output ---------------------\n\n")
            f.write("number of cells in averaging:\n")
            f.write("ncells = %d\n" % len(cells))
            f.write("cumulative area of cells (m^2):\n")
            f.write("area = %.6e\n" % area(cells))

        # 3a. Compute requested integrated quantities (if required)
        # Grab any special fluxes.
        special_fns = {}
        for i in cfg['integrated_outputs']:
            if isinstance(i, tuple):
                special_fns[i[0]] = i[1]
        fluxes = compute_fluxes(cells, cfg['variable_map'], cfg['species'], gmodel, special_fns)
            
        if cfg['output_format'] == 'verbose':
            if len(cfg['integrated_outputs']) > 0:
                print "onedval: Writing out integrated quantities"
                f.write("\n---------------------\n")
                f.write("Integrated quantities\n")
                f.write("---------------------\n")
                f.write("\n")

                for flux in cfg['integrated_outputs']:
                    if isinstance(flux, str):
                        if flux == 'mass flux':
                            f.write("mass flux (kg/s)\n")
                            f.write("m_dot = %.6e\n\n" % fluxes['mass'])
                        elif flux == 'momentum flux':
                            f.write("momentum flux (kg.m/s^2)\n")
                            f.write("mom_dot = %s\n\n" % fluxes['mom'])
                        elif flux == 'energy flux':
                            f.write("energy flux (W)\n")
                            f.write("e_dot = %.6e\n\n" % fluxes['energy'])
                        elif flux == 'species mass flux':
                            for sp in cfg['species']:
                                isp = gmodel.get_isp_from_species_name(sp)
                                f.write("mass flux of %s (kg/s)\n" % sp)
                                f.write("m%s_dot = %.6e\n\n" % (sp, fluxes['species'][isp]))
                        else:
                            print "Requested integrated quantity: ", flux
                            print "is not part of the list of available integrated quantities."
                            print "Bailing out!"
                            sys.exit(1)
                    else:
                        f.write("flux of %s\n")
                        f.write("flux = %.6e\n\n" % fluxes[flux[0]])

        if cfg['output_format'] == 'as_data_file':
            for flux in cfg['integrated_outputs']:
                if isinstance(flux, str):
                    if flux == 'mass flux':
                        int_quants['mass flux'] = fluxes['mass']
                    elif flux == 'momentum flux':
                        int_quants['momentum flux'] = abs(fluxes['mom'])
                    elif flux == 'energy flux':
                        int_quants['energy flux'] = fluxes['energy']
                    elif flux == 'species mass flux':
                        int_quants['species mass flux'] = {}
                        for sp in cfg['species']:
                            isp = gmodel.get_isp_from_species_name(sp)
                            int_quants['species mass flux'][sp] = fluxes['species'][isp]
                    else:
                        print "Requested integrated quantity: ", flux
                        print "is not part of the list of available integrated quantities."
                        print "Bailing out!"
                        sys.exit(1)
                else:
                    int_quants[flux[0]] = fluxes[flux[0]]
        
        # 3b. Compute requested one_d_properties
        if cfg['output_format'] == 'verbose':
            if len(cfg['one_d_averages']) > 0:
                print "onedval: Writing out one-dimensionalised quantities"
                f.write("\n------------------------------\n")
                f.write("One-dimensionalised quantities\n")
                f.write("------------------------------\n")

        phis_all = {}
        for avg in cfg['one_d_averages']:
            if avg == 'area-weighted':
                phis = area_weighted_avg(cells, cfg['one_d_outputs'], cfg['variable_map'])
                phis_all[avg] = copy(phis)
                if cfg['output_format'] == 'verbose':
                    f.write("-- area-weighted average --\n\n")
                    pretty_print_props(f, phis, cfg['species'], cfg['one_d_outputs'])
                    f.write("\n")
            elif avg == 'mass-flux-weighted':
                phis = mass_flux_weighted_avg(cells, cfg['one_d_outputs'], cfg['variable_map'])
                phis_all[avg] = copy(phis)
                if cfg['output_format'] == 'verbose':
                    f.write("-- mass-flux-weighted average --\n\n")
                    pretty_print_props(f, phis, cfg['species'], cfg['one_d_outputs'])
                    f.write("\n")
            elif avg == 'flux-conserved':
                phis, result = stream_thrust_avg(cells, cfg['one_d_outputs'], cfg['variable_map'], cfg['species'], gmodel)
                if result != 'success':
                    print "WARNING: Something went wrong trying to compute flux-conserved averages for slice: ", slice_file
                    if cfg['skip_bad_slices']:
                        print "Skipping this slice and continuing."
                    else:
                        print "Bailing out at this point because 'skip_bad_cells' is set to false."
                        sys.exit(1)

                phis_all[avg] = copy(phis)
                if cfg['output_format'] == 'verbose' and result == 'success':
                    f.write("-- flux-conserved average --\n\n")
                    pretty_print_props(f, phis, cfg['species'], cfg['one_d_outputs'])
                    f.write("\n")
            else:
                print "Requested one-D averaging method: ", avg
                print "is not known or not implemented."
                print "Bailing out!"
                sys.exit(1)

        if cfg['output_format'] == 'as_data_file' and result == 'success':
            A = area(cells)
            pos = avg_pos(cells, cfg['variable_map'])
            data_file_row(f, pos, A, phis_all, int_quants, cfg['one_d_averages'], cfg['one_d_outputs'],
                              cfg['integrated_outputs'], cfg['species'])

    f.close()
    print "onedval: Done."

if __name__ == '__main__':
    main()
    

    

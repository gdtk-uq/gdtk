"""
sensitivity_input_utils.py -- Function to check the sensitivity input 
    dictionary that is pulled in from a config file. Default values 
    are added where appropriate.

.. Author: Luke Doherty (luke.doherty@eng.ox.ac.uk)
           Oxford Thermofluids Institute
           The University of Oxford

  Version: 04-May-2020
"""

import sys, os, copy
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from read_outfile import read_outfile
from perturb_utils import read_case_summary

#--------------------------------------------------------------------------
def sensitivity_input_checker(cfg):
    """Takes the config dictionary and checks that all the necessary
    inputs are in there. Will bail out if important variables are missing,
    and set default values for those that are not crucial.

    Returns the config dictionary after it has been checked over and will
    tell sensitivity to bail out if it finds an issue.
    
    """

    print "Checking sensitivity inputs."

    # Make empty input uncertainties dictionary that we will populate
    #cfg['inputUncertainties'] = {}

    if 'levels' not in cfg:
        cfg['levels'] = 3
        print "The number of perturbation levels to use is not specified."
        print "    Setting it to default value of {0}".format(cfg['levels'])

    if 'outputVariables' not in cfg:
        print "You have not specified a list of which variables are to be analysed."
        print "    Bailing out."
        cfg['bad_input'] = True

    if 'inputVariables' not in cfg:
        print "You have not specified a list of the input variables (that have been perturbed."
        print "    Bailing out."
        cfg['bad_input'] = True

    if 'inputVariablesUncertainty' not in cfg:
        print "You have not specified a dictionary of the input variables uncertainties."
        print "    Bailing out."
        cfg['bad_input'] = True

    if 'typeOfUncertainty' not in cfg:
        print "You have not specified the type of uncertainty (absolute or relative)."
        print "    Bailing out."
        cfg['bad_input'] = True
    
    # Check that we have uncertainties for each input variable
    for variable in cfg['inputVariables']:
        #print variable
        #print cfg['inputVariablesUncertainty'].keys()
        
        if variable not in cfg['inputVariablesUncertainty'].keys():
            print "You have not specified an uncertainty for {0}".format(variable)
            print "    Bailing out."
            cfg['bad_input'] = True
    
    # Read the sensitivity_case_summary file to get the perturbed
    # variables and their various values
    perturbedVariables, dictOfCases = read_case_summary()
    
    # Define the name of the nominal case and load the exit plane data
    nominal = 'case000'
    nominalData, dontNeed = read_outfile('./'+nominal+'/'+cfg['resultsFile'])
    
    # Check for discrepancies between the perturbedVariables list and the
    # inputVariables list
    for var in perturbedVariables:
        if var not in cfg['inputVariables']:
            print "'{0}' was perturbed but does not appear in the 'inputVariables' list.".format(var)
            # Bailing out could be made unnecessary but it would require 
            # signifiant changes to how the main sensitivity function 
            # works since the case number is specific to the variable 
            # position in the inputVariables list. 
            print "    Bailing out."
            cfg['bad_input'] = True
    
    # Clean up the input uncertainty. Not sure if this is better here
    # or in the main function
    if cfg['typeOfUncertainty'] in ['absolute']:
        cfg['inputVariablesUncertainty_abs'] = copy.deepcopy(cfg['inputVariablesUncertainty'])
        cfg['inputVariablesUncertainty_rel'] = {}
        for var in cfg['inputVariables']:
            cfg['inputVariablesUncertainty_rel'][var] = \
                cfg['inputVariablesUncertainty_abs'][var]/\
                dictOfCases['case000'][var]
    else:
        cfg['inputVariablesUncertainty_rel'] = copy.deepcopy(cfg['inputVariablesUncertainty'])
        cfg['inputVariablesUncertainty_abs'] = {}
        for var in cfg['inputVariables']:
            cfg['inputVariablesUncertainty_abs'][var] = \
                cfg['inputVariablesUncertainty_rel'][var]*\
                dictOfCases['case000'][var]
    cfg.pop('inputVariablesUncertainty')
    cfg.pop('typeOfUncertainty')
    
    
    print "Done checking sensitivity inputs."

    return cfg


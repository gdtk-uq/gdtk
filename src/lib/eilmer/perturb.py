#!/usr/bin/env python
""" perturb.py

This script coordinates the running of a series of calculations 
that are perturbations around a nominal condition.

.. Author: Luke Doherty (luke.doherty@eng.ox.ac.uk)
           Oxford Thermofluids Institute
           The University of Oxford
"""

VERSION_STRING = "04-May-2020"

import shlex, subprocess, string
from subprocess import PIPE
import sys, os, gzip
import optparse
from numpy import array, mean, logical_and, zeros
from perturb_input_utils import perturb_input_checker
from perturb_utils import run_command, \
                          set_case_running, set_perturbed_values,\
                          write_case_config, write_case_summary
import copy
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#---------------------------------------------------------------
def run_perturb(cfg):
    """
    Function that accepts a config dictionary and runs perturb.
    """
    #
    cfg['bad_input'] = False
    #
    # First check our input and assign any defaults values
    cfg = perturb_input_checker(cfg)
    #bail out here if there is an issue
    if cfg['bad_input']:
        return -2
    #
    # Unpack the input dictionary of nominal values into the main 
    # configuration dictionary
    for var in cfg['nominalValues'].keys():
        cfg[var] = cfg['nominalValues'][var]
    #
    # Re-assign the dictionary of perturbed values (for later convenience)
    perturbedDict = cfg['perturbedDict']
    #
    if not cfg['createRSA']: # Perturbing for a sensitivity calculation
        # Calculate Nominal condition
        caseString = 'case'+"{0:02}".format(0)+"{0:01}".format(0)
        textString = "Nominal Condition"
        caseDict = copy.copy(cfg)
        caseDict['caseName'] = caseString
        write_case_config(caseDict)
        # Run the nominal case and write the values of the perturbed 
        # variables to a summary file
        set_case_running(caseString, caseDict, textString,'')
        write_case_summary(cfg['variablesToPerturb'],caseDict,caseString,1)
        #
        # Now run all the perturbed conditions
        for k in range(len(cfg['variablesToPerturb'])):
            var = cfg['variablesToPerturb'][k]
            perturbCount = cfg['levels']

            for kk in range(perturbCount):
                if kk != 0:
                    caseString = 'case'+"{0:02}".format(k)+\
                                 "{0:01}".format(kk)
                    textString = var+" perturbed to "+\
                                 str(perturbedDict[var][kk])
                    caseDict = copy.copy(cfg)
                    caseDict[var] = perturbedDict[var][kk]
                    caseDict['caseName'] = caseString
                    # Run the current case
                    set_case_running(caseString, caseDict, \
                                     textString, var)
                    # Write current case to the summary file
                    write_case_summary(cfg['variablesToPerturb'],\
                                       caseDict,caseString,0)
    #
    else: # Perturbing to create a Response Surface
        # Currently we can only fit a response surface through
        # two variables, so we take the first two in the input list
        var1 = cfg['variablesToPerturb'][0] # 'Vs'
        var2 = cfg['variablesToPerturb'][1] # 'pe'

        if cfg['levels'] == 2.5:
            casesToRun = [(2,1),      (1,1),
                                (0,0),
                          (2,2),      (1,2)]
        elif cfg['levels'] == 3:
            casesToRun = [(2,1),(0,1),(1,1),
                          (2,0),(0,0),(1,0),
                          (2,2),(0,2),(1,2)]
        elif cfg['levels'] == 5:
            casesToRun = [(4,3),      (0,3),      (3,3),
                                (2,1),      (1,1),
                          (4,0),      (0,0),      (3,0),
                                (2,2),      (1,2),
                          (4,4),      (0,4),      (3,4)]
            #casesToRun = [      (2,3),      (1,3),
            #              (4,1),      (0,1),      (3,1),
            #                    (2,0),      (1,0),
            #              (4,2),      (0,2),      (3,2),
            #                    (2,4),      (1,4)      ]
        #
        # Run the nominal case first
        caseString = 'case'+"{0:01}{0:01}".format(0,0)
        caseDict = copy.copy(paramDict)
        caseDict['caseName'] = caseString
        textString = "Nominal Case: "+var1+"="+str(perturbedDict[var1][0])+\
                                 "; "+var2+"="+str(perturbedDict[var2][0])
        write_case_config(caseDict)
        set_case_running(caseString, caseDict, textString, '')
        write_case_summary(cfg['variablesToPerturb'], caseDict, caseString, 1)
        #
        # Now run all other cases
        for case in casesToRun:
            if case != (0,0):
                caseString = 'case'+"{0:01}{1:01}".format(case[0],case[1])
                textString = var1+" perturbed to "+\
                             str(perturbedDict[var1][case[0]])+\
                        "\n"+var2+" perturbed to "+\
                             str(perturbedDict[var2][case[1]])
                caseDict = copy.copy(paramDict)
                caseDict['caseName'] = caseString
                caseDict[var1] = perturbedDict[var1][case[0]]
                caseDict[var2] = perturbedDict[var2][case[1]]
                #
                set_case_running(caseString, caseDict, textString, '')
                write_case_summary(cfg['variablesToPerturb'], \
                                   caseDict, caseString, 0)
    #
    return


def main():
    """
    Examine the command-line options to decide the what to do
    and then coordinate a series of calculations using inputs
    that are perturbed around specified nominal values.
    """
    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('-c', '--config_file', dest='config_file',
                  help=("filename for the config file"))
    opt, args = op.parse_args()
    config_file = opt.config_file
    #
    cfg = {}
    #   
    if not cfg: #if the configuration dictionary has not been filled up already, load it from a file
        try: #from Rowan's onedval program
            execfile(config_file, globals(), cfg)
        except IOError as e:
            print "Error {0}".format(str(e))
            print "There was a problem reading the config file: '{0}'".format(config_file)
            print "Check that it conforms to Python syntax."
            print "Bailing out!"
            sys.exit(1)
    #
    run_perturb(cfg)
    #
    return

#---------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "perturb:\n run simulations that are perturbations about a nominal"
        print "   Version:", VERSION_STRING
        print "   To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)

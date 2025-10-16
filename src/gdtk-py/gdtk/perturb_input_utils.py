"""
perturb_input_utils.py -- Function to check the perturb input dictionary
    that is pulled in from a config file. Default values are added 
    where appropriate.

.. Author: Luke Doherty (luke.doherty@eng.ox.ac.uk)
           Oxford Thermofluids Institute
           The University of Oxford

  Version: 04-May-2020
"""

from perturb_utils import set_perturbed_values

def perturb_input_checker(cfg):
    """Takes the config dictionary and checks that all the necessary
    inputs are in there. Will bail out if important variables are missing, 
    and set default values for those that are not crucial.
       
    Returns the config dictionary after it has been checked over and will
    tell perturb to bail out if it finds an issue.
    
    """

    print "Checking perturb inputs."
    #
    if 'runCMD' not in cfg:
        cfg['runCMD'] = './'
        print "Run command not specified. We will assume you are running"
        print "    the simulations on a local machine that does not have"
        print "    a que. Setting to a default run command ('{0}')".format(cfg['runCMD'])
    #
    if 'levels' not in cfg:
        cfg['levels'] = 3
        print "The number of perturbation levels to use is not specified."
        print "    Setting it to default value of {0}".format(cfg['levels'])
    #
    if 'createRSA' not in cfg:
        print "Response surface approximation variable not set."
        print "    Setting it to False."
    #
    if cfg['levels'] == '3-reduced':
        if cfg['createRSA']:
            cfg['levels'] = 2.5
        else:
            print "Setting levels to '3-reduced' is only valid when creating a"
            print "    response surface. Changing this to levels = 3."
            cfg['levels'] = 3
    else:
        cfg['levels'] = int(cfg['levels'])
    #
    if 'variablesToPerturb' not in cfg:
        print "You have not specified a list of which variables are to be perturbed."
        print "    Bailing out."
        cfg['bad_input'] = True
    #
    if 'nominalValues' not in cfg:
        print "You have not specified a dictionary of default values for each variable"
        print "    that will be perturbed."
        print "    Bailing out."
        cfg['bad_input'] = True
    else:
        for var in cfg['variablesToPerturb']:
            if var not in cfg['nominalValues']:
                print "{0} is not present in the nominalValues dictionary."
                print "    Bailing out."
                cfg['bad_input'] = True
    # 
    if 'typeOfPerturbation' not in cfg and 'perturbationMagnitudes' not in cfg:
        print "Neither the typeOfPerturbation or perturbationMagnitudes dictionaries"
        print "    have been specified. Setting a default 'relative' perturbation of"
        print "    1.0 percent for each variable."
        cfg['typeOfPerturbation'] = {}
        cfg['perturbationMagnitudes'] = {}
        for var in cfg['variablesToPerturb']:
            cfg['typeOfPerturbation'][var] = 'relative'
            cfg['perturbationMagnitudes'][var] = [1.0] # This must a list...
    elif 'typeOfPerturbation' not in cfg and 'perturbationMagnitudes' in cfg:
        print "You must specify either both or neither of the typeOfPerturbation"
        print "    and perturbationMagnitudes dictionaries."
        print "    Bailing out."
        cfg['bad_input'] = True
    elif 'typeOfPerturbation' in cfg and 'perturbationMagnitudes' not in cfg:
        print "You must specify either both or neither of the typeOfPerturbation"
        print "    and perturbationMagnitudes dictionaries."
        print "    Bailing out."
        cfg['bad_input'] = True
    else:
        for var in cfg['variablesToPerturb']:
            if var not in cfg['typeOfPerturbation']:
                print "{0} is not present in typeOfPerturbation dictionary.".format(var)
                print "    Bailing out."
                cfg['bad_input'] = True
            if var not in cfg['perturbationMagnitudes']:
                print "{0} is not present in perturbationMagnitudes dictionary.".format(var)
                print "    Bailing out."
                cfg['bad_input'] = True
    # 
    #if 'inputFileTemplate' not in cfg:
    #    print "You have not specified the inputFileTemplate file name."
    #    print "    Bailing out."
    #    cfg['bad_input'] = True
    #
    if 'runScriptTemplate' not in cfg:
        print "You have not specified the runScriptTemplate file name."
        print "    Bailing out."
        cfg['bad_input'] = True
    #
    # Warn the user that only the first two elements of the input 
    # variablesToPerturb list will be perturbed when we want to 
    # create a response surface
    if cfg['createRSA']:
        print "Create Response Surface (RS) option selected."
        if len(cfg['variablesToPerturb'])>2:
            print "NOTE: We can currently only fit a response surface to"
            print "two variables."

    if 'extraFilesToCopy' not in cfg:
        print "No extra files to copy have specified."
        cfg['extraFilesToCopy'] = []

    # Now go through and calculate the perturbed values for each variable
    # of interest.
    if not cfg['bad_input']:
        cfg['perturbedDict'] = {} #this will get populated below
        #
        for variable in cfg['variablesToPerturb']: 
            cfg['perturbedDict'][variable] = set_perturbed_values( \
                                                 cfg['nominalValues'][variable], \
                                                 cfg['perturbationMagnitudes'][variable],\
                                                 cfg['typeOfPerturbation'][variable], \
                                                 cfg['levels'])

            print "Done with perturbed variable {0}.".format(variable)

    print "Done checking perturb inputs."

    return cfg


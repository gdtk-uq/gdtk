"""
Templete Jobfile for use with e4batch.py .

The below is a template jobfile to be used by e4batch.py.
Please erefer to the definition of the global_config dictionary in e4batch.py to
see all the setting options.

Author: Ingo Jahn, Kieran Mackle
Last Modified: 1/9/2020
"""

import numpy as np

# global_config is used to set how the batch will be run.
global_config = {
    """Queue Settings"""
    "QUEUE_TYPE"    : None,  # Allows queue type to be selected for the batch run.
                             # None - will run jobs sequentially on local.
                             # sbatch -  will run via slurm (e.g. for UQ Goliath Cluster)
    "NTASKS"        : 8,     # When using sbatch, sets number of parallel tasks.
    "VERBOSITY"     : 0,     # 0, 1, 2 - set reporting level.

    """File Settings"""
    "BASE_CASE_DIR" : '../simulationsStaticWingFlap/R-A0B0-inv-transient',  # Path to directory
                                        # with reference simulation. This directory will be cloned for
                                        # each case and the modified using parameters listed below.
    "BASE_JOB_FILE" : 'sim_static.lua',  # job filename that needs to be executed to run respective simulation.
    "WORKING_DIR"   : '../TEST/',       # directory where simualtions will be performed and stored
    "CASE_LIST_FILE": 'Case_list.txt',  # File generated during the "prep" stage that lists all the runs
                                        # that are part of the current batch job.
                                        # Can be edited between prep and run stage.
    "CASE_DATA_FILE": 'out_',   # For each run an output file will be generated. These will take
                                # a name like "CASE_DATA_FILECase0.txt", "CASE_DATA_FILECase1.txt", ...
    "DELETE_CASE_DIR": True,    # Flag to define is Case directory will be deleted after a run is completed.
                                # (Usually required to conserve harddrive space)

    """Parametric Variation Settinngs"""
    "VAR_0_NAME"    : 'ALPHA',  # variable VAR_0 name defined in jobfile
    "VAR_1_NAME"    : 'BETA',  # variable VAR_0 name defined in jobfile
    "VAR_REFINE_NAME": 'REFINE',  # varibale name used to set Refinement
    "VAR_REFINE"    : 0.5,  # Level of Refinement for GRID

    "OUTPUT_FILE"   : 'out.txt',  # where simualtion data is stored,
    "OUTPUT_UNITS"  : 'degrees',  # [radians, degrees] sets units of output file for VAR_0 and VAR_1
    }

# Code to generate x-array that defines the parameteric variations.
alpha_list = np.linspace(-1, 1, num=5)
beta_list = np.linspace(-20, 20, num=9)
x = []
for a in alpha_list:
    for b in beta_list:
        x.append([a, b])
x = np.radians(np.array(x))  # convert to radians

# positon dictionary defines parameters that will be varied in the batch run.
position = {
    "VAR_0" : x[:, 0],  # [rad] values for Var_0
    "VAR_1" : x[:, 1],  # [rad] values for Var_1
}

#! /usr/bin/env python3
"""
Code to automatically generate and run batches of Eilmer simulations.

The purpose of this utility is to automate the batch running of Eilmer
simualtions, for example when trying to create aerodynamics tables or performing
parameteric variation studies.
For an example of the jobfile required to run a simulation, please refer to
e4batch-job.py as an example and template.

Three operational modes
1) --prep
    Creates WORKING_DIR and CASE_LIST_FILE. Case list can be edited at this
    stage.
2) --run
    Works through Case list and submits these for processing.
3) --post
    Postprocess and extract data.
4) --prep-run
    Combination of (1) and (2)

Authors: Ingo Jahn, Kieran Mackle
Last Modified: 1 Sept 2020
"""

from getopt import getopt
import numpy as np
import os
import shutil
import sys


class GConf:
    """Define all global-level config vars (with default value if applicable)."""

    _config = {
        "QUEUE_TYPE"    : None,  # If None, will run jobs sequentially.
                                 # If sbatch, will run via slurm (e.g. for UQ Goliath Cluster)
        "NTASKS"        : 8,     # When using sbatch, sets number of parallel tasks.
        "VERBOSITY"     : 0,     # 0, 1, 2 - set reporting level.

        "WORKING_DIR"   : None,  # directory where simualtions will be performed and stored
        "CASE_LIST_FILE": None,  # list containing Cases that will be simulated.    
        "CASE_DATA_FILE": None,  # Name root for output data
        "BASE_CASE_DIR" : None,  # Path to directory for Base simulation
        "BASE_JOB_FILE" : None,  # job filename for Reference simualtion
        "VAR_0_NAME"    : None,  # variable VAR_0 name defined in jobfile
        "VAR_1_NAME"    : None,  # variable VAR_0 name defined in jobfile
        "VAR_REFINE_NAME": None,  # varibale name used to set Refinement
        "VAR_REFINE"    : 1.0,  # Level of Refinement for GRID
        "OUTPUT_FILE"   : None,  # where simualtion data is stored,
        "OUTPUT_UNITS"  : 'radians',  # [radians, degrees] sets units of output file for VAR_0 and VAR_1
        "DELETE_CASE_DIR": False,  # set to True an Case directories will be deleted once run complete.
        }

    class classproperty(property):
        """Class for checking that correct inputs are provided in GConf."""

        def __get__(self, cls, owner):
            return classmethod(self.fget).__get__(None, owner)()

    for var in _config.keys():
        exec("@classproperty\ndef {0}(cls): return cls._config['{0}']".format(var))

    @staticmethod
    def read_inputs(settings):
        """Check if the input variables are valid then overwrit."""
        for option, value in settings.items():
            if option not in GConf._config:
                raise CustomError(
                        "Invalid configuration variable '" + option +
                        "' supplied. Valid configuration variables are [" +
                        ", ".join([GConf._config]) + "]"
                        )
            GConf._config[option] = value

    @staticmethod
    def check_inputs():
        """Checl valid value ranges for each config option."""
        # Config variables not in this dict will NOT be checked
        valid_input_dict = {
            "QUEUE_TYPE"        : [None, 'sbatch'],
            "OUTPUT_UNITS"      : ["radians", "degrees"],
            "DELETE_CASE_DIR"   : [True, False],
            "VERBOSITY"         : [0, 1, 2, 3],
            }
        for option, valid_inputs in valid_input_dict.items():
            if GConf._config[option] not in valid_inputs:
                raise CustomError(
                        "Configuration variable " + option + " = " +
                        str(GConf._config[option]) + " is not valid." +
                        "Valid values are " + str([valid_inputs]))


def createCaseDict_and_update_jobFile(CaseDirName, working_dir_abs, base_dir_abs, GConf, var_0, var_1):
    """
    Create Case specific directory, and update the job file to have correct variables.

    Inputs:
        CaseDirName -
    Output:
        CasePath_abs - absolute path to created Case Directory
    """
    CasePath_abs = os.path.join(working_dir_abs, CaseDirName)
    # delete potentially existing case directory and create new template directory
    if os.path.exists(CasePath_abs):
        shutil.rmtree(CasePath_abs)
    shutil.copytree(base_dir_abs, CasePath_abs)
    print("Created new Case Directory '{}'".format(CasePath_abs))

    # modify the job-file to set VAR_0 and VAR_1
    with open(os.path.join(CasePath_abs, 'temp_job_file.tmp'), 'w+') as new_file:  # create temporary file
        with open(os.path.join(CasePath_abs, GConf.BASE_JOB_FILE), 'r') as old_file:
            for line in old_file:
                if GConf.VAR_0_NAME in line:
                    line = "{0} = {1:8e}\n".format(GConf.VAR_0_NAME, var_0)
                if GConf.VAR_1_NAME in line:
                    line = "{0} = {1:8e}\n".format(GConf.VAR_1_NAME, var_1)
                if GConf.VAR_REFINE_NAME is not None:
                    if GConf.VAR_REFINE_NAME in line:
                        line = "{0} = {1:8e}\n".format(GConf.VAR_REFINE_NAME, GConf.VAR_REFINE)
                new_file.write(line)
    os.remove(os.path.join(CasePath_abs, GConf.BASE_JOB_FILE))
    shutil.move(os.path.join(CasePath_abs, 'temp_job_file.tmp'),
                os.path.join(CasePath_abs, GConf.BASE_JOB_FILE))  # move temp file to replace jobfile

    return CasePath_abs


def main(uo_dict):
    """Run the main piece of code."""

    # Read the jobfile
    exec(open(uo_dict["--job"]).read(), globals())
    GConf.read_inputs(global_config)
    # GConf.check_inputs()

    # set current directory
    start_dir = os.getcwd()

    # get working directory
    if os.sep in GConf.WORKING_DIR:
        relative_dir = GConf.WORKING_DIR.split(os.sep)
        working_dir_abs = os.path.abspath(os.path.join(*relative_dir))
    else:
        relative_dir = '.'
        working_dir_abs = os.path.dirname(os.path.realpath(__file__))

    if "--prep" in uo_dict or "--prep-run" in uo_dict or "--run" in uo_dict:
        # check that BASE_CASE_DIR exists
        if os.sep in GConf.BASE_CASE_DIR:
            relative_dir = GConf.BASE_CASE_DIR.split(os.sep)
            base_dir_abs = os.path.abspath(os.path.join(*relative_dir))
        else:
            relative_dir = '.'
            base_dir_abs = os.path.dirname(os.path.realpath(__file__))
        if not os.path.isdir(base_dir_abs):
            raise CustomError("BASE_CASE_DIR='{}' doesn't exist.".format(base_dir_abs))

    # do preparation
    if "--prep" in uo_dict or "--prep-run" in uo_dict:
        # check if WORKING_DIR exists and/or create
        if os.path.isdir(working_dir_abs):
            print("WORKING_DIR='{}' already exists.".format(working_dir_abs))
        else:
            print("Creating new WORKING_DIR='{}'".format(working_dir_abs))
            os.mkdir(working_dir_abs)

        # check that BASE_JOB_FILE exists
        if not os.path.isfile(os.path.join(base_dir_abs, GConf.BASE_JOB_FILE)):
            raise CustomError("BASE_JOB_FILE='{}' doesn't exist.".format(GConf.BASE_JOB_FILE))

        # check that VAR_NAMEs only occur once in jobfile
        file = open(os.path.join(base_dir_abs, GConf.BASE_JOB_FILE), 'r').read()
        count_0 = file.count(GConf.VAR_0_NAME)
        count_1 = file.count(GConf.VAR_1_NAME)
        if not (count_0 == 1 and count_1 == 1):
            raise CustomError("Check BASE_JOB_FILE. "
                              + "{0} occurs {1:d}; ".format(GConf.VAR_0_NAME, count_0)
                              + "{0} occurs {1:d}; ".format(GConf.VAR_1_NAME, count_1)
                              + "both must occur once only.")
        if GConf.VAR_REFINE_NAME is not None:
            count_refine = file.count(GConf.VAR_REFINE_NAME)
            if not count_refine == 1:
                raise CustomError("Check BASE_JOB_FILE. "
                                  + "{0} occurs {1:d}; ".format(GConf.VAR_REFINE_NAME, count_refine)
                                  + "must occur once only.")

        # create Case_list
        with open(GConf.CASE_LIST_FILE, 'w+') as f:
            f.write("# Case File for bulk_run_CFD.py\n")
            f.write("# jobfile = {}\n".format(uo_dict["--job"]))
            f.write("# working_dir_abs = {}\n".format(working_dir_abs))
            f.write("# base_dir_abs = {}\n".format(base_dir_abs))
            f.write("# base_job_file = {}\n".format(GConf.BASE_JOB_FILE))
            f.write("# pos0:Case pos1:{0} pos2:{1}\n".format(GConf.VAR_0_NAME, GConf.VAR_1_NAME))
            for i, (temp0, temp1) in enumerate(zip(position["VAR_0"], position["VAR_1"])):
                f.write("{0:d} {1:.8e} {2:.8e}\n".format(i, temp0, temp1))
        print("Created CASE_LIST_FILE, '{0}' with {1:d} cases.".format(GConf.CASE_LIST_FILE, i+1))

    # proceed to run simualtions
    if "--prep-run" in uo_dict or "--run" in uo_dict:
        print("")
        if GConf.QUEUE_TYPE is None:
            print("Starting to run jobs sequentially on local machine")
            for i, (var_0, var_1) in enumerate(zip(position["VAR_0"], position["VAR_1"])):  # TODO: replace by option to locae Case_list.txt
                print("")
                print("+++++++++++++++++++++++++++++")
                print("STARTING CASE {}".format(i))
                print("{0}={1}; {2}={3}".format(GConf.VAR_0_NAME, var_0, GConf.VAR_1_NAME, var_1))
                # check if corresponding output file exists
                os.chdir(working_dir_abs)
                outputFileName = GConf.CASE_DATA_FILE+"Case{:04d}.txt".format(i)
                if os.path.isfile(outputFileName):
                    print("Output data file '{}' already exists, skipping ahead.".format(outputFileName))
                else:
                    CaseDirName = "Case{:04d}".format(i)
                    CasePath_abs = createCaseDict_and_update_jobFile(CaseDirName, working_dir_abs, base_dir_abs, GConf, var_0, var_1)
                    os.chdir(CasePath_abs)
                    # print("Changed or current Case Directory='{}'".format(CasePath_abs))

                    # prep:
                    print("")
                    prep_gas_model_string = "prep-gas ../../constants/ideal-air.inp ideal-air-gas-model.lua > LOG-prep-gas"
                    print("START: prep gas model:", prep_gas_model_string)
                    os.system(prep_gas_model_string)
                    print("    DONE: prep gas model.")

                    print("")
                    prep_sim_string = "e4shared --job={0} --prep > LOG-prep".format(GConf.BASE_JOB_FILE.split('.')[0])
                    print("START: prep sim:", prep_sim_string)
                    os.system(prep_sim_string)
                    print("    DONE prep sim.")

                    # run
                    print("")
                    run_sim_string = "e4shared --job={0} --run > LOG-run".format(GConf.BASE_JOB_FILE.split('.')[0])
                    print("START: run sim:", run_sim_string)
                    os.system(run_sim_string)
                    print("    DONE run sim.")

                    # extract load data
                    print("")
                    print("START: run extract load:")
                    tags = ["wing", "flap"]
                    tag_files = []
                    os.chdir(start_dir)
                    for t, tag in enumerate(tags):
                        tag_files.append("{0}-{1}-loads.txt".format(CaseDirName, tag))
                        run_string = ('python3 compute_forces_and_moments.py --job={0} --absolute-path={1} --tindx-plot=all --compute-loads-on-group={2} --output-file={3}'
                                      .format(GConf.BASE_JOB_FILE.split('.')[0], CasePath_abs, tag, tag_files[t]))
                        print("    tag={} ->".format(tag), run_string)
                        os.system(run_string)
                        shutil.move(os.path.join(CasePath_abs,
                                                 "{}".format(tag_files[t])),
                                    os.path.join(working_dir_abs,
                                                 "{}".format(tag_files[t])))
                    print("    DONE extract loads")

                    # write date into Case specific output file
                    os.chdir(working_dir_abs)
                    with open(outputFileName, 'w+') as fnew:
                        fnew.write("# Case{:04d}\n".format(i))
                        fnew.write("# VAR_0:{0}={1}\n".format(GConf.VAR_0_NAME, var_0))
                        fnew.write("# VAR_1:{0}={1}\n".format(GConf.VAR_1_NAME, var_1))
                        fnew.write("# pos0:tag[-] pos1:Time[s] pos2:Fxi[N] pos3:Fyi[N] pos4:Mzi[Nm] pos5:Fxv[N] pos6:Fyv[N] pos7:Mzv[Nm]\n")
                        for file, tag in zip(tag_files, tags):
                            with open(file, 'r') as f:
                                data = f.readlines()
                                fnew.write("{0} {1}".format(tag, data[-1]))  # \n not required as already included in data[-1]

                    if GConf.DELETE_CASE_DIR is True:  # remove the completed Case directories
                        shutil.rmtree(CasePath_abs)
                        print("Deleted Case Directory: {}".format(CasePath_abs))

        elif GConf.QUEUE_TYPE == 'sbatch':
            print("Starting to submit batch jobs using 'sbatch' on slurm queue")
            for i, (var_0, var_1) in enumerate(zip(position["VAR_0"], position["VAR_1"])):  # TODO: replace by option to locae Case_list.txt
                print("")
                print("+++++++++++++++++++++++++++++")
                print("STARTING CASE {}".format(i))
                print("{0}={1}; {2}={3}".format(GConf.VAR_0_NAME, var_0, GConf.VAR_1_NAME, var_1))
                CaseDirName = "Case{:04d}".format(i)

                # check if corresponding output tag-files exist
                tags = ["wing", "flap"]
                tag_files = []
                os.chdir(working_dir_abs)
                tag_exist_flag = 1
                for t, tag in enumerate(tags):
                    tag_files.append("{0}-{1}-loads.txt".format(CaseDirName, tag))
                    if os.path.isfile(tag_files[t]):
                        tag_exist_flag *= 1
                    else:
                        tag_exist_flag *= 0

                if tag_exist_flag == 1:
                    print("Tag files for current case already exists, skipping ahead.")
                else:
                    # create Case Directory and go there.
                    CasePath_abs = createCaseDict_and_update_jobFile(CaseDirName, working_dir_abs, base_dir_abs, GConf, var_0, var_1)
                    os.chdir(CasePath_abs)

                    # create slurm file.
                    slurm_file_name = "prep-run-Case{:04d}.slurm".format(i)
                    with open(slurm_file_name, 'w+') as f:
                        # slurm file header
                        f.write("#!/bin/bash\n")
                        f.write("#SBATCH --job-name={}\n".format(CaseDirName))
                        f.write("#SBATCH --nodes=1\n")
                        f.write("#SBATCH --ntasks={:d}\n".format(GConf.NTASKS))
                        f.write("\n")
                        f.write("module load mpi/openmpi-x86_64\n")
                        f.write("\n")

                        # slurm code to run simualtion
                        prep_gas_model_string = "prep-gas ../../constants/ideal-air.inp ideal-air-gas-model.lua > LOG-prep-gas"
                        f.write(prep_gas_model_string + "\n")
                        prep_sim_string = "e4shared --job={0} --prep > LOG-prep".format(GConf.BASE_JOB_FILE.split('.')[0])
                        f.write(prep_sim_string + "\n")
                        run_sim_string = "e4shared --job={0} --run > LOG-run".format(GConf.BASE_JOB_FILE.split('.')[0])
                        f.write(run_sim_string + "\n")
                        f.write("\n")

                        # slurm code to extract loads and
                        f.write("cd {}\n".format(start_dir))
                        for t, tag in enumerate(tags):
                            run_string = ('python3 compute_forces_and_moments.py --job={0} --absolute-path={1} --tindx-plot=all --compute-loads-on-group={2} --output-file={3}'
                                          .format(GConf.BASE_JOB_FILE.split('.')[0], CasePath_abs, tag, tag_files[t]))
                            f.write(run_string + "\n")
                            f.write("mv {0} {1} \n".format(os.path.join(CasePath_abs,
                                                           "{}".format(tag_files[t])),
                                                           os.path.join(working_dir_abs,
                                                           "{}".format(tag_files[t]))))
                            f.write("\n")
                        # delete Case directory
                        f.write("cd {}\n".format(working_dir_abs))
                        f.write("rm -r {}\n".format(CaseDirName))  # Note empty CaseDir will remain as slurm script executed from within

                    print("slurm file '{0}' has been generated for: {1}".format(slurm_file_name, CaseDirName))
                    os.system("sbatch {}".format(slurm_file_name))
                    print("Job 'sbatch {}' has been submitted".format(slurm_file_name))

        else:
            raise CustomError("Option for GConf.QUEUE_TYPE not supported.")

    if "--post" in uo_dict:  # perform post-processing tasks
        # post-processing routine that combines all the force data into a single file.
        print("")
        print("Starting Postprocessing and colating load data.")

        # write date into Case specific output file
        os.chdir(working_dir_abs)
        outputFileName = GConf.OUTPUT_FILE
        print("Data will be written to: '{}'".format(outputFileName))
        with open(outputFileName, 'w+') as f:
            # write data file header
            f.write("# Force-Moment Lookup tables\n")
            f.write("# BaseCase: {}\n".format(GConf.BASE_CASE_DIR))
            f.write("# BaseJob: {}\n".format(GConf.BASE_JOB_FILE))
            # GConf.ANGLE = 'radians'
            GConf.ANGLE = 'degrees'
            if GConf.ANGLE == 'radians':
                f.write("# pos0:ALPHA[rad] pos1:ALPHA_DOT[rad/s] pos2:BETA[rad] pos3:BETA_DOT[rad/s]"
                        + " pos4:Fx_wing[N] pos5:Fy_wing[N] pos6:Mz_wing[Nm]"
                        + " pos7:Fx_elev[N] pos8:Fy_elev[N] pos9:Mz_elev[Nm]\n")
            elif GConf.ANGLE == 'degrees':
                f.write("# pos0:ALPHA[deg] pos1:ALPHA_DOT[deg/s] pos2:BETA[deg] pos3:BETA_DOT[deg/s]"
                        + " pos4:Fx_wing[N] pos5:Fy_wing[N] pos6:Mz_wing[Nm]"
                        + " pos7:Fx_elev[N] pos8:Fy_elev[N] pos9:Mz_elev[Nm]\n")

            for i, (var_0, var_1) in enumerate(zip(position["VAR_0"], position["VAR_1"])):  # TODO: replace by option to locae Case_list.txt
                CaseDirName = "Case{:04d}".format(i)
                print("")
                print("Processing: {}".format(CaseDirName))
                print("{0}={1}; {2}={3}".format(GConf.VAR_0_NAME, var_0, GConf.VAR_1_NAME, var_1))

                # check if corresponding output tag-files exist
                tags = ["wing", "flap"]
                tag_files = []
                os.chdir(working_dir_abs)
                files_exist = True
                for t, tag in enumerate(tags):
                    tag_files.append("{0}-{1}-loads.txt".format(CaseDirName, tag))
                    if not os.path.isfile(tag_files[t]):
                        print("WARNING: data file = '{}' missing, skipping Case".format(tag_files[t]))
                        files_exist = False
                        break
                if files_exist is False:
                    continue  # skip ahead to look at next Case.

                # create string for current case
                if GConf.ANGLE == 'radians':
                    string = "{} 0.0 {} 0.0".format(var_0, var_1)
                elif GConf.ANGLE == 'degrees':
                    string = "{} 0.0 {} 0.0".format(np.degrees(var_0), np.degrees(var_1))
                for file, tag in zip(tag_files, tags):
                    with open(file, 'r') as ftag:
                        line = ftag.readlines()
                    data = line[-1].rstrip().split(" ")
                    data = [float(d) for d in data]
                    string = string + " {0:.8e} {1:.8e} {2:.8e}".format(data[1]+data[4],
                                                                        data[2]+data[5],
                                                                        data[3]+data[6])
                f.write(string+"\n")
                print("Adding following string to file:\n{}".format(string))

    return 0


class CustomError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def check_inputs(uo_dict):
    """Check all mandatory options have been provided."""
    reqd_options = ["--job"]
    for op in reqd_options:
        if op not in uo_dict:
            raise CustomError("".join(("bulk_run_CFD.py requires argument '",
                              op, "', but this argument was not provided.\n")))

    # Check that jobfile exists
    if not os.path.isfile(uo_dict["--job"]):
        raise CustomError("".join(("Jobfile '", uo_dict["--job"], "' not found,",
                          " check that file path is correct.\n")))


def print_usage():
    print("")
    print("To be completed.")
    print("")


short_options = ""
long_options = ["help", "job=", "prep", "prep-run", "run", "post"]


if __name__ == '__main__':
    user_options = getopt(sys.argv[1:], short_options, long_options)
    uo_dict = dict(user_options[0])

    if len(user_options[0]) == 0 or "--help" in uo_dict:
        print_usage()
        sys.exit(1)

    else:
        check_inputs(uo_dict)
        main(uo_dict)
        print("\n")
        print("SUCCESS.")
        print("\n")

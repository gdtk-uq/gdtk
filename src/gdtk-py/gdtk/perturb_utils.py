"""
perturb_utils.py -- Small untility functions needed by the
                    perturb and sensitivity  programs.

.. Author: Luke Doherty (luke.doherty@eng.ox.ac.uk)
           Oxford Thermofluids Institute
           The University of Oxford

  Version: 04-May-2020
"""

import sys, os
import shlex, subprocess, string
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#-------------------------------------------------------------------------
def run_command(cmdText):
    """
    Run the command as a subprocess. Copied verbatim from nenzfr_utils.py
    """
    
    # Flush before using subprocess to ensure output is in the right order.
    sys.stdout.flush()

    if (type(cmdText) is list):
        args = cmdText
    else:
        args = shlex.split(cmdText)
    print "About to run cmd:", string.join(args)
    # p = subprocess.Popen(args)
    # wait until the subprocess is finished
    # stdoutData, stderrData = p.communicate()
    return subprocess.check_call(args)


def set_case_running(caseString, caseDict, textString, var):
    """
    A short function to set a given case running in an appropriately
    named sub-directory.
    """
    print 60*"-"
    print caseString
    print textString
    #
    if var not in ['coreRfraction','CoreRadiusFraction']:
        # Create sub-directory for the current case
        run_command('mkdir ./'+caseString)
        #
        # Set up the run script
        scriptFileName = caseDict['runScriptTemplate']
        create_file_from_template(scriptFileName, './'+caseString+\
                                  '/'+scriptFileName, caseDict)
        #
        # Set up the input file
        if 'inputFileTemplate' in caseDict:
            inputFileName = caseDict['inputFileTemplate']
            create_file_from_template(inputFileName,'./'+caseString+\
                                      '/'+inputFileName, caseDict)
        #
        # If required, copy additional files to the sub-directory
        for file in caseDict['extraFilesToCopy']:
            command_text = 'cp ./'+file+' ./'+caseString+'/'
            run_command(command_text)
        #
        # Change into the sub-directory, ensure the run script and 
        # input file are exectuable, then execute the run script
        os.chdir(caseString)
        run_command('chmod u+x '+scriptFileName)
        if 'inputFileTemplate' in caseDict:
            run_command('chmod u+x '+inputFileName)
        print ""
        print caseDict['runCMD']+scriptFileName
        print ""
        # I am not sure how to replace the following line with the 
        # run_command function
        os.system(caseDict['runCMD']+scriptFileName)
        os.chdir('../')
    else:
        # Below is very specific to nenzfr. I don't know how best to remove
        # this so I've left it. LukeD 05-05-2020
        if not os.path.exists(caseString):
            # Copy nominal case data into a new case directory
            command_text = 'cp -r case000 '+caseString
            run_command(command_text)
            #
            # Change into new directory and call nenzfr --just-stats
            os.chdir(caseString)
            command_text = \
                'nenzfr.py --just-stats --CoreRadiusFraction='+\
                str(caseDict[var])
                #str(DictOfCases[caseString][perturbedVariables.index(var)])
            print ""
            print command_text
            run_command(command_text)
            print ""
            os.chdir('../')
    #
    return

def create_file_from_template(templateFileName, outPutFileName, substituteDict):
    """Function that will write out a new file by undertaking a dictionary 
    substitution on a nominated template file.
    """
    #
    fp = open(templateFileName, 'r')
    text = fp.read()
    fp.close()
    template = string.Template(text)
    text = template.safe_substitute(substituteDict)
    fp = open(outPutFileName, 'w')
    fp.write(text)
    fp.close()
    #
    return

def set_perturbed_values(NominalValue, PerturbationMagnitude, TypeOfPerturbation, Levels):
    """
    
    """
    if TypeOfPerturbation == "relative":
        delta = [k/100.0*NominalValue for k in PerturbationMagnitude]
    else:
        delta = PerturbationMagnitude
    #
    if Levels == 3:
        if len(delta) == 1:
            # have [delta]
            values = [NominalValue, NominalValue+delta[0],\
                                    NominalValue-delta[0]]
        elif len(delta) == 2:
            # have [+delta, -delta]
            values = [NominalValue, NominalValue+delta[0],\
                                    NominalValue-delta[1]]
    elif Levels == 5:
        if len(delta) == 1:
            # have [delta]
            values = [NominalValue, NominalValue+delta[0],\
                                    NominalValue-delta[0],\
                                    NominalValue+2.0*delta[0],\
                                    NominalValue-2.0*delta[0]]
        elif len(delta) == 2:
            # have [delta, "2"*delta]
            values = [NominalValue, NominalValue+delta[0],\
                                    NominalValue-delta[0],\
                                    NominalValue+delta[1],\
                                    NominalValue-delta[1]]
        elif len(delta) == 3:
            # have [+delta, -delta, "2"*delta]
            values = [NominalValue, NominalValue+delta[0],\
                                    NominalValue-delta[1],\
                                    NominalValue+delta[2],\
                                    NominalValue-delta[2]]
        elif len(data) == 5:
            # have [+delta, -delta, +"2"*delta, -"2"*delta]
            values = [NominalValue, NominalValue+delta[0],\
                                    NominalValue-delta[1],\
                                    NominalValue+delta[2],\
                                    NominalValue-delta[3]]
        
    return values

def write_case_summary(varList,caseDict,caseString,newfile):
    """
    A short function to write the values of the perturbed variables 
    for the current case to a summary file.
    """
    # Define some format strings
    #l = max( max([len(k) for k in varList]), 15 ) 
    #dataFormat = '{0:>'+str(l)+'.6g}'
    #titleFormat = '{0:{fill}>'+str(l)+'}'
    dataFormat = {k : '{0:>'+"{0}".format(max(11,len(k))+2)+".5g}" \
                      for k in varList}
    titleFormat = {k : '{0:{fill}>'+"{0}".format(max(11,len(k))+2)+"}" \
                       for k in varList}
     
    # For the first time, we create a new file and write the header information.
    # Each following case is appended to the existing file.
    if newfile == 1:
        fout = open('perturbation_cases.dat','w')
        # Write title line
        fout.write('{0:>7}'.format('#'))
        for k in tuple(varList):
            fout.write(titleFormat[k].format(k,fill=''))
        fout.write('\n')
        # Underline the title
        for k in varList:
            fout.write(titleFormat[k].format('-',fill='-'))
        fout.write('{0:->7}'.format('-'))
        fout.write('\n')
    else:
        fout = open('perturbation_cases.dat','a')
    # Now write out the data for the current case
    fout.write('{0:>7}'.format(caseString))
    for k in varList:
        fout.write(dataFormat[k].format(caseDict[k]))
    fout.write('\n')
    fout.close()

def read_case_summary(FileToRead=None):
    """
    Reads the file "perturbation_cases.dat" to determine which variables
    have been perturbed and values for the perturbed variables for each
    case consituting the sensitivity calculation.

    :returns: perturbedVariables - a list of variable names
            : DictOfCases - a dictionary with the case names as keys
                            and the values of each perturbed variable.
    """
    if FileToRead is None:
        fp = open('perturbation_cases.dat','r')
    else:
        fp = open(FileToRead,'r')

    varList = fp.readline().strip().split(" ")
    perturbedVariables = [k for k in varList if k!="#" and k!=""]
    fp.readline()
    DictOfCases = {}
    for line in fp.readlines():
        caseData = line.strip().split(" ")
        caseName = caseData[0]
        caseValues = [float(k) for k in caseData if k!=caseName and k!=""]
        DictOfCases[caseName] = dict(zip(perturbedVariables,caseValues))
        #DictOfCases[caseName] = [float(k) for k in caseData \
        #                         if k!=caseName and k!=""]
    fp.close()
    return perturbedVariables, DictOfCases


def write_case_config(caseDict):
    """Short function that writes a file summarising the values of
    an input dictionary. 
    """
    fout = open('nominal_case.config','w')
    for k in range(len(caseDict)):
        fout.write('{0}'.format(caseDict.keys()[k]+':   '))
        fout.write('{0}'.format(caseDict.values()[k]))
        fout.write('\n')
    fout.close()


# The below function was originally used in nenzfr_sensitivity. Parts of
# the code have now been incorporated into the function "set_case_running"
# (see above). LukeD 05-05-2020
# 
#def perturb_CoreRadiusFraction(var, perturbedVariables,\
#                               DictOfCases, levels):
#    """
#    Perturbation of the CoreRadiusFraction may be completed
#    as a post-processing step using the nominal condition
#    flow solution. We don't need to compute separate solutions.
#
#    This function copies the nominal condition solution to a
#    new case file and then runs nenzfr with the --just-stats
#    option and the perturbed CoreRadiusFraction value.
#    """
#
#    for kk in range(levels):
#        if kk != 0:
#            caseString = 'case'+"{0:02}".format(\
#                 perturbedVariables.index(var))+\
#                 "{0:01}".format(kk)
#            # Only if required do we proceed with copying the
#            # nominal case data and calling 'nenzfr.py'
#            if not os.path.exists(caseString):
#                print "Perturbing "+var
#                #command_text = 'rm -r '+caseString
#                #run_command(command_text)
#
#                # Copy nominal case data into a new case
#                # directory
#                command_text = 'cp -r case000 '+caseString
#                run_command(command_text)
#                # Change into new directory and call nenzfr --just-stats
#                os.chdir(caseString)
#                command_text = \
#                    'nenzfr.py --just-stats --CoreRadiusFraction='+\
#                    str(DictOfCases[caseString][perturbedVariables.index(var)])
#                run_command(command_text)
#                os.chdir('../')



#! /usr/bin/env python3
# e4-prep-python.py
"""
Prepare a simulation, writing the block files in batches, concurrently.

Usage:
  $ e4-prep-parallel --job=<jobName> --n-workers=<int> --nb-per-task=<int> \
                     --complex --prep=<prepOption>

Authors:
  PAJ and RJG

Version:
  2020-06-27: original code
  2021-04-21: We now have classic "prep" and new "prep-flow" options.
"""

import sys
import os
from getopt import getopt
import shlex, subprocess
import concurrent.futures


def run_command(cmd):
    """
    Run the command as a subprocess, saving the standard streams to text files.
    Input is a tuple, as indicated below.
    """
    cmdText, tagText = cmd
    if (type(cmdText) is list):
        args = cmdText
    else:
        args = shlex.split(cmdText)
    with open(f"logs/prep-log-stdout-{tagText}.txt","w") as out, \
         open(f"logs/prep-log-stderr-{tagText}.txt","w") as err:
        completedCmd = subprocess.run(args, stdout=out, stderr=err)
    return completedCmd.returncode

def main(jobName, nWorkers, nBlocksPerTask, doComplexPrep=False, prepOption="prep-flow"):
    """
    Do the real work with the user-supplied parameters.
    """
    print(f"jobName={jobName}, nWorkers={nWorkers}, nBlocksPerBatch={nBlocksPerTask}")
    #
    inputScriptName = jobName
    if prepOption == "prep-flow": inputScriptName += "-flow"
    inputScriptName += ".lua"
    print("Input script file name: %s" % inputScriptName)
    if not os.path.exists(inputScriptName):
        raise RuntimeError(f"Input script not found: {inputScriptName}")
    if not os.path.exists("logs"): os.mkdir("logs")
    prepProgName = "e4shared"
    if doComplexPrep: prepProgName = "e4zshared"
    cmdTxt = f"{prepProgName} --job={jobName} --{prepOption} --no-block-files"
    tagTxt = f"{0}"
    print(f"cmd={cmdTxt}, tag={tagTxt}")
    returnCode = run_command((cmdTxt, tagTxt))
    if returnCode != 0:
        raise RuntimeError(f"Preparation of config files failed. cmd={cmdTxt}, tag={tagTxt}")
    #
    global workers
    if nWorkers > 1:
        workers = concurrent.futures.ThreadPoolExecutor(max_workers=nWorkers)
    else:
        workers = None
    listFile = open(f"config/{jobName}.list", 'r')
    lastLine = listFile.readlines()[-1]
    nBlocks = int(lastLine.strip().split()[0]) + 1
    print(f"nBlocks={nBlocks}")
    cmds = []
    for i in range(0, nBlocks, nBlocksPerTask):
        n1 = i
        n2 = min(i+nBlocksPerTask, nBlocks)
        cmdTxt = f"{prepProgName} --job={jobName} --{prepOption} --no-config-files --only-blocks=\"{n1}..<{n2}\""
        tagTxt = f"{n2}"
        print(f"cmdTxt={cmdTxt}, tagTxt={tagTxt}")
        cmds.append((cmdTxt, tagTxt))
    if workers:
        print("Build block files concurrently.")
        my_futures = [workers.submit(run_command, cmd) for cmd in cmds]
        returnCodes = [fut.result() for fut in my_futures]
        workers.shutdown()
    else:
        print("Build block files serially.")
        returnCodes = [run_command(cmd) for cmd in cmds]
    if any(returnCodes):
        raise RuntimeError(f"One or more prep commands failed. returnCodes={returnCodes}")
    return

#--------------------------------------------------------------------

shortOptions = "hf:w:b:cp:"
longOptions = ["help", "job=", "n-workers=", "nb-per-task=", "complex", "prep"]

def printUsage():
    print("Prepare a simulation, writing the block files in batches, concurrently.")
    print("Usage: e4-prep-parallel [--help | -h] [--job=<jobName> | -f <jobName>]\n" +
          "                        [--n-workers=<int> | -w <int>] [--nb-per-task=<int> | -b <int>]\n" +
          "                        [--complex | -c] [--prep=<option> | -p <option>]")
    print("")
    return

if __name__ == '__main__':
    print("Begin e4-prep-parallel...")
    #
    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) == 0 or \
           "--help" in uoDict or \
           "-h" in uoDict:
        printUsage()
    else:
        if "--job" in uoDict:
            jobName = uoDict.get("--job", "")
        elif "-f" in uoDict:
            jobName = uoDict.get("-f", "")
        else:
            raise Exception("Job name is not specified.")
        jobName, ext = os.path.splitext(jobName)
        if "--n-workers" in uoDict:
            n_workers = int(uoDict.get("--n-workers", "1"))
        elif "-w" in uoDict:
            n_workers = int(uoDict.get("-w", "1"))
        else:
            n_workers = 1
        if "--nb-per-task" in uoDict:
            nbpt = int(uoDict.get("--nb-per-task", "1"))
        elif "-b" in uoDict:
            nbpt = int(uoDict.get("-b", "1"))
        else:
            nbpt = 1
        if "--complex" in uoDict:
            doComplexPrep = True
        elif "-c" in uoDict:
            doComplexPrep = True
        else:
            doComplexPrep = False
        if "--prep" in uoDict:
            prepOption = uoDict.get("--prep", "prep-flow")
        elif "-p" in uoDict:
            prepOption = uoDict.get("-p", "prep-flow")
        else:
            prepOption = "prep-flow"
        main(jobName, n_workers, nbpt, doComplexPrep, prepOption)
    print("Done.")
    sys.exit(0)

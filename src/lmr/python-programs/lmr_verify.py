#!/usr/bin/env python3
# A small program to perform a verification test for lmr5.
#
# The verification test involves running the same calculation
# on grids of increasing refinement. This program coordinates
# the preparation, running and results assembly for a 
# sequence of grids of various refinment levels.
#
# Output:
#    Goes in subdirectory case-[<number>|<tag>]
#
# Author: Rowan J. Gollan
# Date: 2023-07-15
#
# History: 2023-08-19
#          Re-work into a program that can be used in any directory

import click
import yaml
import os
import subprocess
import shutil
from math import log

filesToCopyFilename = "FILES_TO_COPY"

@click.command()
@click.option("-cf", "--case-file", "caseFile",
              default="cases.yml",
              help="YAML file with verification cases.",
              show_default=True)
@click.option("-cn", "--case-number", "caseNumber", type=int,
              default=0,
              help="Case number to process.",
              show_default=True)
@click.option("-ct", "--case-tag", "caseTag",
              default="none",
              help="""
\b
Case tag to process.
If supplied, this is used in preference to a supplied case-number.
""",
              show_default=True)
@click.option("--norms", "normsStr",
              default="rho",
              help="Error norms for calculation as comma-separated list.",
              show_default=True)
@click.option("-glf", "--grid-level-file", 'gridLevelFile',
              default="grid-levels.py",
              show_default=True,
              help="""
\b
A python file to define the various levels of grid refinement.
At end of execution of the python file, a list of dictionaries
called 'gridLevels' is expected to be configured.
"""
              )
@click.option("-gl", "--grid-levels", 'levelsToExec',
              default="-1",
              show_default=True,
              help="""
\b
A comma-separated list of grid levels to include in the verification.
For example, to process grid levels k=[0,1,2], use
   --grid-levels="0,1,2"
The value of "-1" (default) indicated to process *all* grid levels
set in the grid levels file.
"""
              )
@click.option("-cd", "--common-dir", "commonDir",
              default="common",
              show_default=True,
              help="""
\b
The name of the directory containing files that are common to all runs
of the verification. Inside this directory a file called FILES_TO_COPY
is expected. This file is plain text with one filename per line.
The filenames listed in FILES_TO_COPY are used in all case directories
for the verification runs.
"""
              )
@click.option("-ppo", "--post-process-only", "postProcessOnly",
              is_flag=True,
              default=False,
              show_default=True,
              help="Only perform post-processing portion (norm and observed order calcs)"
              )
def verify(caseFile, caseNumber, caseTag, normsStr, gridLevelFile, levelsToExec, commonDir, postProcessOnly):
    """Run a verification test.

    \b
    Example for grid levels 2, 1 and 0 with density and temperature norms:
       > lmr-verify --grid-levels="0,1,2" --norms="rho,T"
    Example for norms of density and velocities for case with tag 'ausmdv-2nd-order':
       > lmr-verify --norms="rho,vel.x,vel.y" -ct ausmdv-2nd-order

    """
    # 0. Prepare command-line options
    case = findCase(caseFile, caseNumber, caseTag)
    if case == None:
        print("No case to process found. Exiting.")
        exit(1)
    
    levelsToExec = [int(lvl) for lvl in levelsToExec.split(",")]
    levelsToExec.sort(reverse=True)
    norms = normsStr.split(",")
    norms = [norm.strip() for norm in norms]

    # 1. Prepare gridLevels and simFiles
    exec(open(gridLevelFile).read(), globals())
    if len(levelsToExec) == 1 and levelsToExec[0] == -1:
        levelsToExec = list(reversed(range(len(gridLevels))))

    if not postProcessOnly:
        simFiles = assembleSimFiles(commonDir)
        # 2. Execute verification runs
        prepareGridLevels(case, gridLevels, levelsToExec, commonDir, simFiles)
        runGridLevels(case, levelsToExec)

    # 3. Post-process results
    computeNorms(case, levelsToExec, normsStr)
    assembleResults(case, gridLevels, levelsToExec, norms)

def assembleSimFiles(commonDir):
    with open(commonDir + "/" + filesToCopyFilename, "r") as f:
        simFiles = [line.strip() for line in f.readlines() if line != ""]
    return simFiles


def findCase(caseFile, caseNumber, caseTag):
    if caseTag != "none":
        return findCaseByTag(caseFile, caseTag)

    else:
        return findCaseByNumber(caseFile, caseNumber)

def findCaseByTag(caseFile, caseTag):
    with open(caseFile, "r") as f:
        docs = list(yaml.safe_load_all(f))
        for case in docs:
            if case["tag"] == caseTag:
                case["label"] = caseTag
                return case

    # If we get here, we failed to find case tag
    print(f"Unable to find tag= {caseTag} in file= '{caseFile}'")
    return None

def findCaseByNumber(caseFile, caseNumber):
    with open(caseFile, "r") as f:
        docs = list(yaml.safe_load_all(f))
        if caseNumber < len(docs):
            docs[caseNumber]["label"] = f"{caseNumber:03d}"
            return docs[caseNumber]
        else:
            print(f"Unable to find {caseNumber=} in file= '{caseFile}'")
            print(f"Number of cases found is: {len(docs)} (cases count from 0)")
            return None

def prepareGridLevels(case, gridLevels, levelsToExec, commonDir, simFiles):
    print("Preparing cases.")
    cwd = os.getcwd()
    caseDir = "case-" + case["label"]
    for k in levelsToExec:
        subDir = caseDir + f"/k-{k}"
        print(f"-- Grid level {k=}")
        os.makedirs(subDir, exist_ok=True)
        os.chdir(subDir)
        with open("config.txt", "w") as f:
            f.write(buildConfigStr(case, gridLevels[k]))
        with open("run.sh", "w") as f:
            f.write(buildRunStr())
        for f in simFiles:
            fname = f"../../{commonDir}/{f}"
            shutil.copy(fname, ".")
        os.chdir(cwd)

def buildConfigStr(case, grid):
    cfgStr = ""
    for param, value in case.items():
        cfgStr += f"{param} = {value!r}\n"
    for param, value in grid.items():
        if param == "dx": continue
        cfgStr += f"{param} = {value!r}\n"
    return cfgStr

def buildRunStr():
    return (
        f"lmr prep-grid\n"
        f"lmr prep-flow\n"
        f"lmr run-steady\n"
        )

def runGridLevels(case, levelsToExec):
    print("Running simulations with manufactured solution source terms.")
    cwd = os.getcwd()
    caseDir = "case-" + case["label"]
    for k in levelsToExec:
        print(f"-- Grid level {k=}")
        subDir = caseDir + f"/k-{k}"
        os.chdir(subDir)
        cmd = "sh run.sh"
        proc = subprocess.run(cmd.split(), capture_output=True, text=True)
        if not checkRun(proc):
            raise RuntimeError(f"There was a problem running grid level {k}")
        os.chdir(cwd)
    return

def checkRun(proc):
    lines = proc.stdout.split("\n")
    reason = None
    for line in lines:
        if line.find("STOP-REASON") != -1:
            reason = line.split()[1]
    if reason != "relative-global-residual-target":
        return False
    return True

def computeNorms(case, levelsToExec, normsStr):
    print("Computing norms.")
    cwd = os.getcwd()
    caseDir = "case-" + case["label"]
    for k in levelsToExec:
        print(f"-- Grid level {k=}")
        subDir = caseDir + f"/k-{k}"
        os.chdir(subDir)
        cmd = f"lmr compute-norms -f -n {normsStr} -r ref-soln.lua -o ../error-norms-k-{k}.txt"
        proc = subprocess.run(cmd.split())
        assert proc.returncode == 0, f"Failed to compute norms for grid level {k=}"
        os.chdir(cwd)
    return

def assembleResults(case, gridLevels, levelsToExec, norms):
    print("Assembling output files.")
    L1 = {}; L2 = {}; Linf = {}; dx = {}
    caseDir = "case-" + case["label"]
    f = open(f"{caseDir}/error-norms-{caseDir}.dat", "w")
    header = "# dx "
    for norm in norms:
        header += f"L1-{norm}  L2-{norm}  Linf-{norm} "
    header += "\n"
    f.write(header)
    for k in levelsToExec:
        fname = f"{caseDir}/error-norms-k-{k}.txt"
        with open(fname, "r") as file:
            doc = yaml.safe_load(file)
        L1[k] = {}; L2[k] = {}; Linf[k] = {};
        dx[k] = gridLevels[k]['dx']
        row = f"{dx[k]:20.12e} "
        for norm in norms:
            L1[k][norm] = doc[norm]["L1"]
            L2[k][norm] = doc[norm]["L2"]
            Linf[k][norm] = doc[norm]["Linf"]
            row += f"{L1[k][norm]:20.12e} {L2[k][norm]:20.12e} {Linf[k][norm]:20.12e} "
        row += "\n"
        f.write(row)
    f.close()
    # Observed order calculation
    if len(levelsToExec) == 1:
        # We can't extract an observed order from one grid refinement.
        # Warn the user and exit.
        print(f"Only one refinement level requested with k={levelsToExec[0]}")
        print("So observed order calculation is not available.")
        return
    f = open(f"{caseDir}/observed-order-{caseDir}.dat", "w")
    header = "# dx "
    for norm in norms:
        header += f"L1-{norm}  L2-{norm}  Linf-{norm} "
    header += "\n"
    f.write(header)
    for i in range(len(levelsToExec)-1):
        kp1 = levelsToExec[i]
        k = levelsToExec[i+1]
        logr = log(dx[kp1]/dx[k])
        row = f"{dx[k]:20.12e} "
        for norm in norms:
            pL1 = log(L1[kp1][norm]/L1[k][norm])/logr
            pL2 = log(L2[kp1][norm]/L2[k][norm])/logr
            pLinf = log(Linf[kp1][norm]/Linf[k][norm])/logr
            row += f"{pL1:20.12e} {pL2:20.12e} {pLinf:20.12e} "
        row += "\n"
        f.write(row)
    f.close()
    return

if __name__ == '__main__':
    verify()


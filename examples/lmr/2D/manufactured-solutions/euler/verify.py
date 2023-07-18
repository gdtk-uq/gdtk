#!/usr/bin/env python
# A small program to perform a verification test for lmr5.
#
# The verification test involves running the same calculation
# on grids of increasing refinement. This program coordinates
# the preparation, running and results assembly for a 
# sequence of grids of various refinment levels.
#
# Usage:
# $ python verify.py --norms="rho,vel.x" --number-cells="8,16,32"
#
# Output:
#    Goes in subdirectory case-[<number>|<tag>]
#
# Author: Rowan J. Gollan
# Date: 2023-08-15
#

import click
import yaml
import os
import subprocess
import shutil
from math import log

domain_length = 1.0
simFiles = ['job.lua',
            'constants.txt',
            'very-viscous-air.lua',
            'analytic_solution.py',
            'udf-source-terms.lua',
            'udf-bc.lua',
            'ref-soln.lua',
            'fill-fn.lua']

commonDir = 'common-files'

@click.command()
@click.option("-cf", "--case-file", "caseFile",
              default="verification-cases.yml",
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
@click.option("-nc", "--number-cells", 'number_cells',
              default="8,16,32",
              show_default=True,
              help="""
\b
List of number of cells for the various levels of grid refinement.
A single integer is provided for each grid because discretisation
is assumed to be equal in the x- and y-directions.
The values are sorted smallest to largest and grids are processed
in that sorted order.
"""
              )
def verify(caseFile, caseNumber, caseTag, normsStr, number_cells):
    """Run a verification test.

    \b
    Example for grids 8x8, 16x16, 32x32 with density and temperature norms:
       > python verify.py --norms="rho,T" --number-cells="8,16,32"
    Example for norms of density and velocities on 16x16, 32x32, 64x64:
       > python verify.py --norms="rho,vel.x,vel.y" -nc 16,32,64

    """
    # 0. Prepare command-line options
    case = findCase(caseFile, caseNumber, caseTag)
    if case == None:
        print("No case to process found. Exiting.")
        exit(1)
    
    ncellsList = [int(ncell) for ncell in number_cells.split(",")]
    ncellsList.sort()
    norms = normsStr.split(",")
    norms = [norm.strip() for norm in norms]

    # 1. Execute verification runs
    prepareGridLevels(case, ncellsList)
    runGridLevels(case, ncellsList)

    # 2. Post-process results
    computeNorms(case, ncellsList, normsStr)
    assembleResults(case, ncellsList, norms)



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


def prepareGridLevels(case, ncellsList):
    print("Preparing cases.")
    cwd = os.getcwd()
    ngrids = len(ncellsList)
    caseDir = "case-" + case["label"]
    for i, ncells in enumerate(ncellsList):
        k = ngrids - i
        subDir = caseDir + f"/k-{k}"
        print(f"-- Grid level {k=}: {ncells}x{ncells}")
        os.makedirs(subDir, exist_ok=True)
        os.chdir(subDir)
        with open("config.txt", "w") as f:
            f.write(buildConfigStr(case, ncells))
        with open("run.sh", "w") as f:
            f.write(buildRunStr())
        for f in simFiles:
            fname = f"../../{commonDir}/{f}"
            shutil.copy(fname, ".")
        os.chdir(cwd)


def buildConfigStr(case, ncells):
    cfgStr = ""
    for param in case.keys():
        cfgStr += f"{param} = {case[param]!r}\n"
    cfgStr += f"ncells = {ncells}\n"
    return cfgStr

def buildRunStr():
    return (
        f"lmr prep-grid\n"
        f"lmr prep-flow\n"
        f"lmr run-steady\n"
        )

def runGridLevels(case, ncellsList):
    print("Running simulations with manufactured solution source terms.")
    cwd = os.getcwd()
    ngrids = len(ncellsList)
    caseDir = "case-" + case["label"]
    for i, ncells in enumerate(ncellsList):
        k = ngrids - i
        print(f"-- Grid level {k=}: {ncells}x{ncells}")
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

def computeNorms(case, ncellsList, normsStr):
    print("Computing norms.")
    cwd = os.getcwd()
    ngrids = len(ncellsList)
    caseDir = "case-" + case["label"]
    for i, ncells in enumerate(ncellsList):
        k = ngrids - i
        print(f"-- Grid level {k=}: {ncells}x{ncells}")
        subDir = caseDir + f"/k-{k}"
        os.chdir(subDir)
        cmd = f"lmr compute-norms -f -n {normsStr} -r ref-soln.lua -o ../error-norms-k-{k}.txt"
        proc = subprocess.run(cmd.split())
        assert proc.returncode == 0, f"Failed to compute norms for grid level {k=}"
        os.chdir(cwd)
    return

def assembleResults(case, ncellsList, norms):
    print("Assembling output files.")
    L1 = {}; L2 = {}; Linf = {}; dx = {}
    ngrids = len(ncellsList)
    caseDir = "case-" + case["label"]
    f = open(f"{caseDir}/error-norms-{caseDir}.dat", "w")
    header = "# dx "
    for norm in norms:
        header += f"L1-{norm}  L2-{norm}  Linf-{norm} "
    header += "\n"
    f.write(header)
    for i, ncells in enumerate(ncellsList):
        k = ngrids - i
        fname = f"{caseDir}/error-norms-k-{k}.txt"
        with open(fname, "r") as file:
            doc = yaml.safe_load(file)
        L1[k] = {}; L2[k] = {}; Linf[k] = {};
        dx[k] = domain_length/ncells
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
    if len(ncellsList) == 1:
        # We can't extract an observed order from one grid refinement.
        # Warn the user and exit.
        print(f"Only one refinement level requested with ncells={ncellsList[0]}")
        print("So observed order calculation is not available.")
        return
    f = open(f"{caseDir}/observed-order-{caseDir}.dat", "w")
    header = "# dx "
    for norm in norms:
        header += f"L1-{norm}  L2-{norm}  Linf-{norm} "
    header += "\n"
    f.write(header)
    for k in range(len(ncellsList)-1,0,-1):
        logr = log(dx[k+1]/dx[k])
        row = f"{dx[k]:20.12e} "
        for norm in norms:
            pL1 = log(L1[k+1][norm]/L1[k][norm])/logr
            pL2 = log(L2[k+1][norm]/L2[k][norm])/logr
            pLinf = log(Linf[k+1][norm]/Linf[k][norm])/logr
            row += f"{pL1:20.12e} {pL2:20.12e} {pLinf:20.12e} "
        row += "\n"
        f.write(row)
    f.close()
    return

if __name__ == '__main__':
    verify()


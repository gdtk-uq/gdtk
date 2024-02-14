#!/usr/bin/env python3
"""
Continuous integration script for Eilmer 5.

@author: NNG and RJG
@date: 2024-02-14

History:
    2024-02-14 Code borrowed heavily from Nick's cibot.py as a starting point.
"""

import yaml
import os
import subprocess
import sys
import pytest

tmpIgnoreList = ['--ignore=2D/convex-corner',
                 '--ignore=2D/wedge']

cfgFile = 'autobot-cfg.yaml' 
cfg = {}

def readConfig():
    with open(cfgFile, "r") as f:
        global cfg
        cfg = yaml.safe_load(f)
    return


def logrun(cmd):
    process = subprocess.run(cmd.split(), capture_output=True, text=True)
    sys.stdout.write(cmd + " ---> STDOUT: \n")
    sys.stdout.write(process.stdout)
    sys.stdout.write("\n")
    sys.stdout.write(cmd + " ---> STDERR: \n")
    sys.stdout.write(process.stderr)
    sys.stdout.write("\n")
    sys.stdout.flush()
    return process.returncode, process.stdout

def pullCode():
    os.chdir(os.path.join(os.getenv("HOME"), cfg['lmrDir']))
    # 0. clean
    subprocess.run("make clean".split(), check=True)
    # 1. prepare for pull, declare where we are presently
    oldhead = subprocess.check_output('git rev-parse HEAD'.split()).decode(sys.stdout.encoding)
    print("Pulling code, old head was: ", oldhead)
    # 2. Now get the new code
    cmd = "git pull"
    flag, out = logrun(cmd)
    print("Pull complete, output is:\n", out)
    # 3. find out what we got!?!
    newhead = subprocess.check_output('git rev-parse HEAD'.split()).decode(sys.stdout.encoding)
    isnew = newhead!=oldhead
    revid = subprocess.check_output('git rev-parse --short HEAD'.split()).decode(sys.stdout.encoding)
    if isnew:
        print("Found new commit(s): {}",format(newhead.strip()))
    return isnew, revid

def buildCode(buildCmd):
    os.chdir(cfg['lmrDir'])
    # 0. clean
    subprocess.run("make clean".split(), check=True)
    # 1. now build, Bob!
    flag, out = logrun(buildCmd)
    isgood = flag==0
    return isgood, out

def installCode(bldCmd):
    os.chdir(cfg['lmrDir'])
    subprocess.run((bldCmd + ' install').split(), check=True)
    return

def runTests(revid):
    os.chdir(cfg['testDir'])
    print(os.getcwd())
    # 0. clean up
    process = subprocess.run("rm -rf lmr".split(), check=True)
    # 1. prepare area
    cmd = f"rsync -r {cfg['exDir']} ./"
    subprocess.run(cmd.split(), check=True)
    # 2. Launch
    os.chdir("lmr")
    logFile = cfg['logFilePrefix'] + '-' + revid + '.log'
    args = tmpIgnoreList + [f"--log-file={logFile}"]
    returncode = pytest.main(args)
    return returncode







    

readConfig()
print(cfg)

bldCmd = "make"
#buildCode(bldCmd)

#installCode(bldCmd)



revid = subprocess.check_output('git rev-parse --short HEAD'.split()).decode(sys.stdout.encoding)
print(revid)
runTests(revid.strip())





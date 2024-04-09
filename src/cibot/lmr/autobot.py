#!/usr/bin/env python3
"""
Continuous integration script for Eilmer 5.

@author: NNG and RJG
@date: 2024-02-14

History:
    2024-02-14 Code borrowed heavily from Nick's cibot.py as a starting point.
"""

import click
import yaml
import os
import subprocess
import sys
import pytest
from datetime import datetime
from time import sleep
import shlex

cfgFile = 'autobot-cfg.yaml'
cfg = {}
receiveList = []
logfile = ""
mailfile = ""
DRY_RUN = False

@click.command()
@click.option("--dry-run", "dryRun",
              is_flag=True,
              default=False,
              help="Perform a dry run.",
              show_default=True
)
@click.option("--ignore-new", "ignoreNew",
              is_flag=True,
              default=False,
              help="Ignore whether code is new or not, just launch anyway.",
              show_default=True
)
def launchTests(dryRun, ignoreNew):
    global DRY_RUN
    DRY_RUN = dryRun
    readConfig()
    global mailfile
    mailfile = f"{cfg['logDir']}/{cfg['mailFilePrefix']}-{datetime.now().strftime('%Y-%m-%d+%H-%M-%S')}.txt"
    isnew, revid = pullCode()
    if isnew or ignoreNew:
        buildAndTest(revid)
    return 0

def compileAddressList():
    global DRY_RUN
    addresses = cfg['devsOnReceiveList']
    os.chdir(cfg['lmrDir'])
    committer = subprocess.check_output('git log -1 --pretty=format:%ae'.split(), text=True)
    if committer not in addresses:
        addresses.append(committer)
    if DRY_RUN:
        addresses = [cfg['testEmail']]
    return addresses

def mailout(subject, message):
    global mailfile, DRY_RUN
    with open(mailfile, "w") as f:
        f.write(message)
    # pause a little while, just to give a network
    # file system enough time to get the file written
    sleep(1)
    addresses = compileAddressList()
    recv = ",".join(addresses)
    cmd = f"mail -s '{subject}' {recv} < {mailfile}"
    if DRY_RUN:
        print(cmd)
        return 0
    else:
        # just do it!
        flag = os.system(cmd)
        # pause another few seconds
        # we need the backgrounded mail process
        # to finish before the main process finishes
        sleep(3)
        return flag

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
    os.chdir(cfg['lmrDir'])
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
    revid = revid.strip()
    if isnew:
        print("Found new commit(s): {}",format(newhead.strip()))
    return isnew, revid

def buildCode():
    os.chdir(cfg['lmrDir'])
    # 0. clean
    subprocess.run("make clean".split(), check=True)
    # 1. now build, Bob!
    cmd = cfg['buildCommand']
    process = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
    isgood = process.returncode==0
    return isgood, process.stdout, process.stderr

def whenBuildIsBad(output, error):
    os.chdir(cfg['lmrDir'])
    commithash = subprocess.check_output('git log -1 --pretty=format:%h'.split(), text=True)
    subject = f"lmr-autobot> Build failed in commit {commithash} ({datetime.now()})"

    commitinfo = subprocess.check_output('git log -1'.split(), text=True)
    body = f"lmr-autobot> Build failure detected in the following commit:\n\n"
    body += '\t'.join(('\n'+commitinfo).splitlines(True))
    body += "\n"
    body += "lmr-autobot> STDERR:\n\n"
    body += error
    body += "\n"
    body += "lmr-autobot> STDOUT:\n\n"
    body += output

    return subject, body

def installCode():
    os.chdir(cfg['lmrDir'])
    cmd = f"{cfg['buildCommand']} INSTALL_DIR={cfg['instDir']} install"
    subprocess.run(shlex.split(cmd), check=True)
    return

def runTests(revid):
    os.chdir(cfg['testDir'])
    # 0. clean up
    process = subprocess.run("rm -rf lmr".split(), check=True)
    # 1. prepare area
    cmd = f"rsync -r {cfg['exDir']} ./"
    subprocess.run(cmd.split(), check=True)
    # 2. Launch
    os.chdir("lmr")
    cmd = "pytest"
    process = subprocess.run(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    global logfile
    logfile = f"{cfg['logDir']}/{cfg['logFilePrefix']}-{revid}.log"
    with open(logfile, "w") as f:
        f.write(process.stdout)

    return process.returncode

def checkTests(returncode):
    msg = f"lmr-autobot> pytest returned with exit code: {pytest.ExitCode(returncode)}\n"

    exitMsg = {
        pytest.ExitCode.OK: "Build and tests successful.",
        pytest.ExitCode.TESTS_FAILED: "Failure(s) when running tests.",
        pytest.ExitCode.INTERRUPTED: "pytest was interrupted during testing.",
        pytest.ExitCode.INTERNAL_ERROR: "An internal error prevented pytest from completing.",
        pytest.ExitCode.USAGE_ERROR: "pytest was used incorrectly.",
        pytest.ExitCode.NO_TESTS_COLLECTED: "No tests were collected by pytest."
    }

    msg += f"lmr-autobot> {exitMsg.get(pytest.ExitCode(returncode))}\n"

    os.chdir(cfg['lmrDir'])
    commithash = subprocess.check_output('git log -1 --pretty=format:%h'.split(), text=True)
    commitinfo = subprocess.check_output('git log -1'.split(), text=True)
    if returncode == 0:
        subject = f"lmr-autobot> Build successful; all tests passed for commit {commithash} ({datetime.now()})"
        msg += "lmr-autobot> Integration tests passed for the following commit:\n\n"
    else:
        # all other codes are some form of failure
        subject = f"lmr-autobot> Test failure(s) in commit {commithash} ({datetime.now()})\n"
        msg += f"lmr-autobot> Failure in tests for the following commit:\n\n"
    msg += '\t'.join(('\n'+commitinfo).splitlines(True))
    msg += '\n'
    msg += f"lmr-autobot> Output from pytest is:\n\n"
    with open(logfile, "r") as f:
        msg += f.read()

    return subject, msg


def buildAndTest(revid):
    isgood, out, err = buildCode()
    if isgood:
        installCode()
        returncode = runTests(revid)
        subject, message = checkTests(returncode)
    else:
        subject, message = whenBuildIsBad(out, err)

    flag = mailout(subject, message)
    return 0

if __name__ == '__main__':
    launchTests()





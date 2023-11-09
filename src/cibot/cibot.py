#!/usr/bin/python3
"""
Python code for running continuous integration tests of Eilmer 4.

Notes:

@author: Nick Gibbons

Version history:
   2022-02-15 RJG converted internal references for use on cloudbot with user account testbot
"""

import fcntl
import subprocess
from datetime import datetime
import os
import sys
from re import findall
from mail import *

USER = 'testbot'

SENDER = "eilmer4bot@gmail.com"
HANDLER = 'n.gibbons@uq.edu.au'
DEVS = ['n.gibbons@uq.edu.au',
        'lachlan.whyborn@uqconnect.edu.au',
        'kyle.damm@uqconnect.edu.au',
        'p.jacobs@uq.edu.au',
        'r.gollan@uq.edu.au']
INSTALL_DIR = f'/home/{USER}/gdtkinst'
BUILD = "make DMD=ldc2 WITH_MPI=1 WITH_COMPLEX_NUMBERS=1 WITH_NK=1 WITH_CHT=1 WITH_E4DEBUG=1 FLAVOUR=fast INSTALL_DIR={}".format(INSTALL_DIR)
LOCKFILE = f'/home/{USER}/lockfile.lock'
LMR_SRC_DIR = f'/home/{USER}/source/gdtk/src/eilmer'
TEST_DIR = f'/home/{USER}/tests'
LOG_DIR = f'/home/{USER}/tests/logs'
PASSWD_FILE = f'/home/{USER}/.gmail_password'

class Lock(object):
    """ A basic system-wide lock. Prevents two version of this script from running at once. """
    def __enter__(self):
        self.fp = open(LOCKFILE)
        fcntl.flock(self.fp.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)

    def __exit__(self, _type, value, tb):
        fcntl.flock(self.fp.fileno(), fcntl.LOCK_UN)
        self.fp.close()

class LogFileIO(object):
    """ Redirect STDOUT and STDERR to a file object for hands off logging """
    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        self.oldstdout = sys.stdout
        self.oldstderr = sys.stderr
        self.logfile = open(self.filename, "w")
        sys.stdout = self.logfile
        sys.stderr = self.logfile

    def __exit__(self, _type, value, tb):
        sys.stdout = self.oldstdout
        sys.stderr = self.oldstderr
        self.logfile.close()

def logrun(command):
    """
    Run a shell command and keep a log of the stdout and stderr combined together.

    Args:
        command: The shell command to be run, already split into arguments.
                 eg. command = ["rsync", "-rvP", "file.txt", "file.bak"]
    Returns:
        exitCode: The exitcode of the command, generally zero if all went well.
        output:   The stdout and stderr of the command, as a wall text

    Notes:
        I tried having separate streams for error and stdout but this is really
        complicated because it requires asynchronous io. So this will have to do.
    """
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    output = []
    while True:
        # Note that in this version of python, nextline is a bytes object
        nextline = process.stdout.readline()
        if nextline == b'' and process.poll() is not None:
            break
        nexttext = nextline.decode(sys.stdin.encoding) # convert bytes to string
        sys.stdout.write(nexttext)
        sys.stdout.flush()
        output.append(nexttext)

    o,e = process.communicate() # This blocks until command is done
    exitCode = process.returncode
    output = ''.join(output)

    return exitCode, output

def pull_code():
    os.chdir(LMR_SRC_DIR)
    subprocess.run("make clean".split(), check=True)

    oldhead = subprocess.check_output('git rev-parse HEAD'.split()).decode(sys.stdout.encoding)
    print("Pulling code, old head was: ", oldhead)
    code, out = logrun("git pull".split())
    newhead = subprocess.check_output('git rev-parse HEAD'.split()).decode(sys.stdout.encoding)
    print("Pull complete, output is:\n", out)
    isnew = newhead!=oldhead
    if isnew: 
        print("Found new commit(s): {}".format(newhead.strip()))
    return isnew

def build_code(buildcommand):
    os.chdir(LMR_SRC_DIR)
    subprocess.run("make clean".split(), check=True)
    exitcode, out = logrun(buildcommand.split())
    isgood = exitcode==0
    return isgood, out

def install_code(buildcommand):
    os.chdir(LMR_SRC_DIR)
    subprocess.run((buildcommand+' install').split(), check=True)

def run_tests():
    os.chdir(TEST_DIR)
    subprocess.run("rm -rf eilmer".split(), check=True)
    subprocess.run("rsync -r ../source/gdtk/examples/eilmer ./".split(), check=True)

    os.chdir(f'{TEST_DIR}/eilmer')
    exitcode, out = logrun('./eilmer-test.rb')
    return out

def check_tests(output):
    """
    Take the stdout from the ruby script and parse it to see whether the tests are okay

    Args:
        output: The stdout and stderr from the eilmer-test.rb script.

    Returns:
        failures: A dictionary of test results, containing the names of each bad test
                  as the keys, and the number of failed tests as the values.

    Notes:
    """
    tests = findall('cmd= ruby [\S]+.rb[\s\S]+?assertions/s', output)

    failures = {}
    for test in tests:
        name = findall('cmd= ruby ([\S]+).rb', test)[0]
        ntests = int(findall('([0-9]+) tests,', test)[0])
        nfailures = int(findall('([0-9]+) failures,', test)[0])
        if nfailures!=0:
            failures[name] = nfailures

    return failures

def extract_execution_time(output):
    """
    Take the stdout from the ruby script and find the total execution time inside it 

    Args:
        output: The stdout and stderr from the eilmer-test.rb script.

    Returns:
        execution_time: string representing the test run time in H:M:S

    Notes:
    """
    for line in output.splitlines():
        if line.startswith('Elapsed time:'):
            time = line.split(':')[-1].strip('.').strip()
            return time
    else:
        print("Warning: No time found in execution output!")
    return

def compile_email_addressees():
    addressees = [email for email in DEVS]

    os.chdir(LMR_SRC_DIR)
    commiter = subprocess.check_output('git log -1 --pretty=format:%ae'.split(), text=True)
    if commiter not in DEVS:
        addressees.append(commiter)
    return addressees

def build_is_bad(buildlog):
    """
    If a failed compile is detected, extract the info of the failed commit and email some people.

    Args:
        buildlog: The stdout and stderr from the compile process.

    Notes: This function is not very clever: It assumes that any line with the word "Error" in it
           is stderr. It would be nice if we could capture the two streams separately but alas.
    """
    addressees = compile_email_addressees()

    os.chdir(LMR_SRC_DIR)
    commithash = subprocess.check_output('git log -1 --pretty=format:%h'.split(), text=True)
    subject = "CI build failure in commit {} ({})".format(commithash, datetime.now())
    print(subject)

    commitinfo = subprocess.check_output("git log -1".split(), text=True)
    body = ['Build failure detected in the following commit:\n',
            commitinfo,
            'Reported stderr was:\n']
    body.extend([line for line in buildlog.splitlines() if 'Error' in line])
    body = '\n'.join(body)

    send_mail(addressees, subject, body)

def tests_are_bad(failedtests):
    """
    If some failed tests are detected, extract the info of the failed tests and email some people.

    Args:
        failedtests: A dictionary of failed test results, returned by check_tests

    Notes:
    """
    addressees = compile_email_addressees()

    os.chdir(LMR_SRC_DIR)
    commithash = subprocess.check_output('git log -1 --pretty=format:%h'.split(), text=True)
    subject = "CI test failure in commit {} ({})".format(commithash, datetime.now())
    print(subject)

    commitinfo = subprocess.check_output("git log -1".split(), text=True)
    body = ['Integration test failure detected in the following commit:\n',
            commitinfo,
            'Failed tests are:']
    body.extend(["   {:18s} : {:3d} failures".format(k,v) for k,v in failedtests.items()])
    body = '\n'.join(body)

    send_mail(addressees, subject, body)

def tests_are_good(execution_time=None):
    """
    Some actions to take if the tests are good.

    Notes: One day I might remove this, I just like seeing the script run : )
    """
    addressees = compile_email_addressees()

    os.chdir(LMR_SRC_DIR)
    commithash = subprocess.check_output('git log -1 --pretty=format:%h'.split(), text=True)
    subject = "CI tests successful for commit {} ({})".format(commithash, datetime.now())
    print(subject)

    commitinfo = subprocess.check_output("git log -1".split(), text=True)
    body = ['Integration tests passed for the following commit:\n',
            commitinfo]
    if execution_time != None: body.append('Execution Time: {}'.format(execution_time))
            
    body = '\n'.join(body)

    send_mail(addressees, subject, body)

def send_mail(addressees, subject, body):
    """
    Instantiate the class in sendmail.py that knows how to send emails.

    Args:
        subject: A string containing the subject line for the email
        body:    The text of the email
    Notes:
        This function just prints a failure message of the email can't be sent. Should be an exception?
    """

    print("In newest send_mail method")
    message=TEMPLATE.format(body=body,handler=HANDLER)
    to = ','.join(addressees)
    create_message_and_send(SENDER, to, subject, message)
    print("Complete!")

def wet_build_and_run():
    """
    Attempt to build the code and then call wet_run to run the tests

    Notes:
        The activity is "wet" in the sense of having consequences, as opposed to a "dry run"
    """
    buildgood, output = build_code(BUILD)
    if buildgood:
        install_code(BUILD)
        wet_run()
    else:
        build_is_bad(output)

def wet_run():
    """
    Assuming a successful build go ahead and run the integration tests.

    Notes:
        The activity is "wet" in the sense of having consequences, as opposed to a "dry run"
    """
    output = run_tests()
    failedtests = check_tests(output)
    nfailed = sum(failedtests.values())

    if nfailed==0:
        execution_time = extract_execution_time(output)
        tests_are_good(execution_time)
    else:
        tests_are_bad(failedtests)

def main(ignorenew=False):
    print("Checking for updates...")
    isnew = pull_code()
    if isnew or ignorenew:
        now = datetime.now().strftime("%Y-%m-%d+%H-%M-%S")
        logfile = f'{LOG_DIR}/{now}.log'
        with LogFileIO(logfile): 
            wet_build_and_run()
    else:
        print("Up to date. Exiting.")
    print("Done...")

if __name__=='__main__':
    print("Hello world!")
    ignorenew=False
    print(sys.argv)
    if len(sys.argv)>1:
        if sys.argv[1]=='ignorenew': ignorenew=True
    #with Lock():
    #    print("Lock aquired...")
    main(ignorenew)

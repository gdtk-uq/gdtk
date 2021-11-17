#! /usr/bin/env python3
# eilmer-test-parallel.py
# Smoke-test the Eilmer flow solver in parallel.
# PJ 2021-11-17, adapted from the Ruby test runner
#

import shlex
import subprocess
import time
import datetime
import os
import os.path
import concurrent.futures

result = subprocess.run(['which', 'gpmetis'], capture_output=True)
gpmetis_exe = result.stdout
found_gpmetis = len(gpmetis_exe) > 0

f = open('test-names.txt', 'r')
lines = f.readlines()
f.close
test_scripts = []
for txt in lines:
    name = txt.strip()
    if ((name.find("metis") < 0) or found_gpmetis):
        test_scripts.append(name)

# print(f"test_scripts={test_scripts}")

def run_test_script(ts):
    work_dir = os.path.dirname(ts)
    text = f"{time.asctime(time.localtime())} {os.getcwd()}\n"
    extension = os.path.splitext(ts)[-1]
    if extension == ".rb":
        cmd = "ruby"
    elif extension in [".test", ".tcl"]:
        cmd = "tclsh"
    else:
        text += "Dodgy extension for test script.\n"
        cmd = ""
    if len(cmd) > 0:
        cmd += " " + os.path.basename(ts)
        text += f"cwd={work_dir} cmd={cmd}\n"
        result = subprocess.run(shlex.split(cmd), cwd=work_dir, capture_output=True)
    return text+result.stdout.decode()

time_start = time.time()

nWorkers = 2
if nWorkers > 1:
    workers = concurrent.futures.ThreadPoolExecutor(max_workers=nWorkers)
    my_futures = [workers.submit(run_test_script, ts) for ts in test_scripts]
    returned_texts = [fut.result() for fut in my_futures]
    workers.shutdown()
    for txt in returned_texts: print(txt)
else:
    for ts in test_scripts:
        returned_text = run_test_script(ts)
        print(returned_text)

time_end = time.time()
print(f"{time.asctime(time.localtime())} Finished tests.")
delta = time_end - time_start
h = int(delta // 3600)
m = int((delta - h*60) // 60)
s = int(delta - h*3600 - m*60)
print(f"Elapsed time: {h} hr, {m} min, {s} sec.")

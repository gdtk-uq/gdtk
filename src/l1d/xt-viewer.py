#! /usr/bin/env python3
# xt-viewer.py
# Display xt-data for the l1d flow simulation.
# The data has previously been accumulated into JSON format files,
# one for each gas slug.
#
# Usage:
# xt-data.py [options] *.json
#
# PJ, First code 2021-01-14
#
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Display xt-data for the l1d4.')
parser.add_argument(dest='fileNames', metavar='filename', nargs='*')
parser.add_argument('-v', '--variable', metavar='varName',
                    dest='varName', action='store',
                    choices={'p', 'rho', 'T', 'vel'}, default='p',
                    help='flow variable name')
parser.add_argument('-l', '--take-log', dest='takeLog', action='store_true',
                    help='take log10 of variable')
parser.add_argument('-x', '--x-range', metavar='x0:dx:x1',
                    dest='xRange', action='store',
                    help='min, delta and max for x')
parser.add_argument('-t', '--t-range', metavar='t0:dt:t1',
                    dest='tRange', action='store',
                    help='min, delta and max for t')
parser.add_argument('-r', '--v-range', metavar='v0:dv:v1',
                    dest='vRange', action='store',
                    help='min, delta and max for variable')
args = parser.parse_args()

print('fileNames=', args.fileNames)
print('varName=', args.varName)
print('takeLog=', args.takeLog)
print('xRange=', args.xRange)
print('tRange=', args.tRange)
print('vRange=', args.vRange)

varName = args.varName if args.varName else 'p'

slugs = []
for fname in args.fileNames:
    with open(fname, 'r') as fp:
        content = fp.read()
    slugs.append(json.loads(content))
    print(slugs[-1].keys())

# Get a single Axes object so that we can put several contour patches onto it.
# We need to explicitly set the levels so that they are the same in all patches.
fig1, ax = plt.subplots()
if args.xRange:
    items = args.xRange.split(':')
    x0 = float(items[0])
    dx = float(items[1])
    x1 = float(items[2])
    print('x0=', x0, 'dx=', dx, 'x1=', x1)
    ticklocs = np.arange(x0, x1+dx/10, dx)
    plt.xticks(ticklocs)
    plt.xlim(x0, x1)
if args.tRange:
    items = args.tRange.split(':')
    t0 = float(items[0])
    dt = float(items[1])
    t1 = float(items[2])
    print('t0=', x0, 'dt=', dt, 't1=', t1)
    ticklocs = np.arange(t0, t1+dt/10, dt)
    plt.yticks(ticklocs)
    plt.ylim(t0, t1)
if args.vRange:
    items = args.vRange.split(':')
    v0 = float(items[0])
    dv = float(items[1])
    v1 = float(items[2])
    print('v0=', v0, 'dv=', dv, 'v1=', v1)
    vLevels = np.arange(v0, v1+dv/10, dv)
else:
    vLevels = None

plts = []
for s in slugs:
    x = np.array(s['x'])
    t = np.array(s['t'])
    v = np.array(s[varName])
    plts.append(ax.contourf(x, t, v, levels=vLevels))

# Assume that all slugs have the same list of variables and units and
# just pick the names out of the first one.
allVarNames = slugs[0]["varNames"]
unitNames = slugs[0]["varUnits"]
ax.set_xlabel(f'x, {unitNames[0]}')
ax.set_ylabel(f't, {unitNames[1]}')
varUnits = unitNames[allVarNames.index(varName)]
ax.set_title(f'xt-diagram of {varName} ({varUnits})')
plt.show()

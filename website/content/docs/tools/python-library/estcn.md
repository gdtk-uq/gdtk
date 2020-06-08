---
title: ESTCN
weight: 30
draft: false
---

## Equilibrium Shock Tunnel Conditions, with Nozzle

This program can be used to estimate flow conditions
for shock-processed flows typical of high-performance
shock-tunnels and expansion tubes.

The gas is assumed to remain in thermochemical equilibrium
and the flow processing is done in decoupled quasi-one-dimensional
wave processes such as shock waves and expansion fans.
For the reflected shock tunnel, this means that the initial,
quiescent test gas is first processed by the incident shock
and subsequently by the reflected shock.
The incident shock sets the inflow conditions for the reflected shock
but there is no further interaction.

The program can do a number of calculations:

* flow in a reflected shock tube with or without a nozzle
* pitot pressure from free-stream flow condition
* stagnation (total) condition from free-stream condition
* code surface condition from free-stream condition

### Example shock-tunnel condition calculation

When run as an application, this program takes its input as
command line arguments, performs the requested calculations and outputs
the gas-state results.
Which particular input parameters you need to supply depends on the
chosen task, however, a typical flow condition for the T4 shock tunnel
with the Mach 4 nozzle may be computed using:

    $ estcn --task=stn --gas=cea-lut-air.lua \
            --T1=300 --p1=125.0e3 --Vs=2414 --pe=34.37e6 --ar=27.0

This condition has an enthalpy of 5.43 MJ/kg and the nozzle-exit condition
has a pressure of 93.6 kPa and a static temperature of 1284 degrees K,
with a flow speed of 2.95 km/s.

### Getting help

To see what specific inputs are required, start the program as:

    $ estcn --help


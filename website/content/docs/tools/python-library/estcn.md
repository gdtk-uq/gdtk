---
title: ESTCN
weight: 30
draft: false
---

## Equilibrium Shock Tunnel Conditions, with Nozzle

This program can be used to estimate flow conditions
for shock-processed flows typical of high-performance
shock-tunnels and expansion tubes.

The program can do a number of calculations from the command line:

* flow in a reflected shock tube with or without a nozzle
* pitot pressure from free-stream flow condition
* stagnation (total) condition from free-stream condition
* code surface condition from free-stream condition

## Getting started

`estcn` is a Python3 program built upon the core gas models and the loadable library
that is wrapped around those models.
To install the estcn program and the loadable library,
move to the `gas` source directory and use the make utility.

    cd dgd/src/gas
    make install

This will also install files associated with the gas models.
Note that you should have your Linux environment set up as described on the
[Quick Start Guide]({{< relref "quick-start" >}}).


## Documentation
- [ESTCN Manual (HTML)](/html/estcn-manual.html)


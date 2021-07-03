---
title: NENZF1D
weight: 40
draft: false
---

## Nonequilibrium nozzle flow

This program can be used to estimate flow conditions
for reflected shock-tunnels when the finite-rate chemical effects
in the nozzle expansion are expected to be important.


## Getting started

`nenzf1d` is a D-language program built upon the core gas-model and flow modules
of our gas-dynamics tool kit.
To install the `nenzf1d` program, move to its source directory
and use the `make` utility.

    cd dgd/src/nenzf1d
    make install

You will also need to install files associated with the gas models.

    cd dgd/src/gas
    make install

Note that you should have your Linux environment set up as described on the
[Quick Start Guide]({{< relref "quick-start" >}}) and have the NASA-Lewis CEA2
program compiled and installed where the gas-model module can find it.


## Documentation
- [NENZF1D Manual (HTML)](/html/nenzf1d-manual.html)



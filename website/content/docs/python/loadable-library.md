---
title: "Loadable library for Python"
description: "ESTCN"
lead: ""
date: 2021-05-26
lastmod: 2021-05-26
draft: false
images: []
menu:
  docs:
    parent: "python"
weight: 10
toc: true
---

The loadable library comes with `eilmer` as a top-level namespace.
Within that namespace, there are packages for gas models,
flow analysis, geometry construction and general numerical methods.
There are also a number of application programs that may be used
from the command line.

## Gas-dynamic functions

This package provides access to functions for the thermochemical gas model
and, built on top of that, simple state-to-state and
stream-tube flow analysis functions.

## Geometry functions
The collection of geometric functions is built to look like the Lua functions
for constructing geometric objects in Eilmer.

## Application programs
Programs built upon the loadable library:
+ [ESTCN]({{< relref "estcn-manual-for-hugo" >}}),
  Calculator for equilibrium shock tunnel conditions, with nozzle.
+ [NENZF1d]({{< relref "nenzf1d-manual-for-hugo" >}}),
  Calculation of Nonequilibrium Nozzle Flow with finite-rate chemistry.

## Getting started

To use the library, it first needs to be built and installed from within
the source code directory for the gas models.

    cd dgd/src/gas
    make install

Note that this command will also install the Python application programs
to the `${DGD}/bin` directory.

If you have not already done so,
you should set the following environment variables
so that the Python interpreter can find the library.

    export PYTHONPATH=${PYTHONPATH}:${DGD}/lib
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DGD}/lib


Remember to refresh your current shell (or log out and log back in) so
that your newly configured environment is available.

## Documentation
+ [Gas-dynamic Package Reference Manual (HTML)]({{< relref "gd-manual-for-hugo" >}})
+ [Geometry Package Reference Manual (HTML)]({{< relref "geometry-manual-for-hugo" >}})
+ [ESTCN, Equilibrium Shock Tunnel Condition, with Nozzle]({{< relref "estcn-manual-for-hugo">}})
+ [NENZF1d, Shock Tunnel Condition, with Nonequilibrium Nozzle]({{< relref  "nenzf1d-manual-for-hugo">}})
+ [Numerical-methods package]({{< relref "nm-manual-for-hugo" >}})







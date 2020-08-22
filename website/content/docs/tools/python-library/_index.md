---
title: Python library and programs
weight: 40
---
## Loadable library for Python
The loadable library comes with `eilmer` as top-level namespace.
Within that namespace, there are packages for gas models,
flow analysis, geometry construction and general numerical methods.
There are also a number of application programs that may be used
from the command line.

### Gas-dynamic functions

This package provides access to functions for the thermochemical gas model
and, built on top of that, simple state-to-state and
stream-tube flow analysis functions.

### Geometry functions
The collection of geometric functions is built to look like the Lua functions
for constructing geometric objects in Eilmer.

### Application programs
Programs built upon the loadable library:
+ [IMOC]({{< relref "moc.md" >}}),
  an interactive program for the isentropic method for characteristics
  for 2-dimensional ideal gas flows.
+ [ESTCN]({{< relref "estcn.md" >}}),
  Calculator for equilibrium shock tunnel conditions, with nozzle.
+ [build-uniform-lut]({{< relref "build-lut.md" >}}),
  Assemble a look-up-table gas model, based on CEA2 calculations
  for a gas mixture in thermochemical equilibrium.
+ [NENZF2]({{< relref "nenzf2.md" >}}),
  Calculation of Nonequilibrium Nozzle Flow with finite-rate chemistry.
+ [Poshax]({{< relref "poshax.md" >}}),
  Calculation of Post-Shock Relaxing flow.

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

- [Gas-dynamic Package Reference Manual (HTML)](/html/library-reference-manual.html)
- [Geometry Package Reference Manual (HTML)](/html/geometry-reference-manual.html)
- [ESTCN, Equilibrium Shock Tunnel Condition, with Nozzle](/html/estcn-manual.html)
- [build-uniform-lut, for building look-up-table gas models](/html/build-lut-manual.html)
- [Numerical-methods package](/html/nm-reference-manual.html)


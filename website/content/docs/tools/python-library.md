---
title: Python library
weight: 40
---
## Loadable library for Python
The loadable library comes with `eilmer` as top-level namespace.
Within that namespace, there are packages for gas models,
flow analysis, geometry construction and general numerical methods.

### Gas-dynamic functions

This package provides access to functions for the thermochemical gas model
and, built on top of that, simple state-to-state and
stream-tube flow analysis functions.

### Geometry functions
The collection of geometric functions is built to look like the Lua functions
for constructing geometric objects in Eilmer.


## Getting started

To use the library, it first needs to be built and installed from within
the source code directory for the gas models.

    cd dgd/src/gas
    make install

You should set the following environment variables, also.

    export PYTHONPATH=${PYTHONPATH}:${DGD}/lib
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DGD}/lib


Remember to refresh your current shell (or log out and log back in) so
that your newly configured environment is available.


## Documentation

- [Gas-dynamic Package Reference Manual (HTML)](/html/library-reference-manual.html)
- [Geometry Package Reference Manual (HTML)](/html/geometry-reference-manual.html)



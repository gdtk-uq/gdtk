---
title: Python library
weight: 40
---

# Gas-dynamic tool-kit for Python

This library provides access to functions for the thermochemical gas model
and, built on top of that, simple state-to-state and
stream-tube flow analysis functions.
There is also a collection of geometric functions.

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

- [Reference Manual (HTML)](/html/library-reference-manual.html)



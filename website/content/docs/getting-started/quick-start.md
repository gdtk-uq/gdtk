---
title: "Quick Start Guide"
date: 2020-09-18
weight: 10
draft: false
---
# Quick Start Guide

This introduction takes you through to running a simple simulation
of supersonic flow over a sharp cone.
It should take less than half an hour to do the full set up and run.


## Prerequisites: background

This section is about you, the user of the Gas Dynamics Toolkit.
We assume that your mathematics, science or engineering background
adequately prepares you for computational fluid dynamics (CFD) analysis.
In particular, we assume that you have a working knowledge of geometry, calculus, mechanics,
and thermo-fluid-dynamics, at least to a second- or third-year university level.
With the Toolkit code, we try to make the analysis of compressible, reacting flow accessible
and reliable; we cannot make it trivial.


## Prerequisites: software

Our main development environment is Linux but the programs can be deployed on
Linux, flavours of Unix such as MacOS-X, and MS-Windows using WSL2.

The core Eilmer, and L1d solvers and their modules are mainly written in the
D programming language for speed and the benefits of compile-time checking.
The pre- and post-processing modes make use of the Lua scripting language
so that we get flexibility and convenient customization.
There is also some Ruby a little Tcl/Tk used in the automated testing scripts.

To run simulations, you need an executable versions of the Eilmer and/or L1d programs.
You may build these executable programs from the source code, as described below.


## Prerequisites: building from source

The main requirement is a D language compiler.
We recommend using the latest stable release of the LLVM D compiler.

To build Eilmer and other programs in the toolkit, you will require:

  + D compilers
      + Binary releases for the latest stable release of the LLVM D compiler (`ldc2` and `ldmd2`)
        may be found at: https://github.com/ldc-developers/ldc/releases
      + An install guide for the LDC compiler is available [here]({{< relref installing-ldc >}}).
  + A C compiler
      + GNU compiler is a good option and comes standard on most systems.
  + The gfortran compiler (and 32-bit libraries)
      + gfortran and gfortran-multilib on Debian/Ubuntu/Mint
      + gcc-gfortran on RedHat/CentOS/Fedora
  + git (to clone the repository)
  + readline development package:
      + libreadline-dev on Debian/Ubuntu/Mint
      + readline-devel on RedHat/CentOS/Fedora
  + ncurses development package:
      + libncurses5-dev on Debian/Ubuntu/Mint
      + ncurses-devel on RedHat/CentOS/Fedora
  + openmpi development package:
      + libopenmpi-dev on Debian/Ubuntu/Mint
      + openmpi-devel on RedHat/CentOS/Fedora
        (after install on RedHat-family systems, load with `module load mpi/openmpi-x86_64`,
        and you might like to place that in your `.bashrc` file so that it's loaded every
        time you start a session)
  + plotutils development package:
      + libplot-dev on Debian/Ubuntu/Mint
      + plotutils-devel on RedHat/CentOS/Fedora (for CentOS 8.x, enable PowerTools repo)
  + foreign-function interface packages for Python and Ruby:
      + python3-cffi on Debian/Ubuntu/Mint and RedHat/CentOS/Fedora
      + ruby-ffi on Debian/Ubuntu/Mint

Additionally, if you want to run the test suite, you will require:

  + Ruby package
  + TCL package
  + the Python sympy package

For viewing and plotting results, we recommend:

  + Paraview
  + Gnuplot

The source code of the Lua interpreter is included in the source code repository.

## Getting the source code

The full source code for the toolkit programs, including a set of examples,
can be found in a public repository on [bitbucket](https://bitbucket.org/cfcfd/dgd-git).
To get your own copy, use the git revision control client to clone the repository
with something like the following command:

    git clone https://bitbucket.org/cfcfd/dgd-git dgd

and within a couple of minutes, depending on the speed of your network connection,
you should have your own copy of the full source tree and the complete repository history.


## Installing Eilmer

The default installation directory is `$HOME/dgdinst`.
To compile and install Eilmer, move into the eilmer source
area and use `make` to coordinate the compiling and installing:

    cd dgd/src/eilmer
    make install

If you are on a Mac, you'll need to give the `make` command an
extra hint:

    make PLATFORM=macosx install


## Installing the loadable library (optional)

The loadable library consists of the gas module functions compiled into
a dynamically loadable library and a collection of Python and Ruby modules
for loading that library.
Once installed, you should be able to use the gas functions from within
a Python or Ruby interpreter.

To compile and install the library and its supporting wrappers,
move to the gas source directory and use `make` again.

    cd dgd/src/gas
    make install

Note that the loadable library needs to be built with the DMD64 compiler
and that you need the Foreign-Function-Interface extensions for your
Python and Ruby interpreters.
On a LinuxMint system these packages are `python-cffi` and `ruby-ffi`.


## Installing L1d (optional)

L1d installs into the same location as Eilmer but, its source code
is in a different location.
To compile and install L1d, move into its source area and, again,
use `make` to coordinate the compiling and installing:

    cd dgd/src/l1d
    make install


## Environment variables

We'll assume you are happy using the default install area `$HOME/dgdinst`.
The next step is to configure your environment to use Eilmer.
You will need to set the
`DGD` variable to point to the top of the installation tree, and
the `DGD_REPO` variable to point to the top of the repository tree.
Note that the installation tree and repository tree are separate. You then
also need to set $PATH, $DGD_LUA_PATH and $DGD_LUA_CPATH to point to
the appropriate places. Some example lines from a .bashrc file are:

    export DGD=$HOME/dgdinst
    export DGD_REPO=$HOME/dgd
    export PATH=$PATH:$DGD/bin
    export DGD_LUA_PATH=$DGD/lib/?.lua
    export DGD_LUA_CPATH=$DGD/lib/?.so


To use the gas models as a loadable library in Python and Ruby,
you should set the following environment variables, also.

    export PYTHONPATH=${PYTHONPATH}:${DGD}/lib
    export RUBYPATH=${RUBYPATH}:${DGD}/lib
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DGD}/lib


Remember to refresh your current shell (or log out and log back in) so
that your newly configured environment is available.


## Running your first Eilmer simulation

To test that everything has worked, you can exercise the flow
solver to simulate the supersonic flow over a 20-deg cone.

    cd ~
    cd dgd/examples/eilmer/2D/sharp-cone-20-degrees/sg
    prep-gas ideal-air.inp ideal-air-gas-model.lua
    e4shared --prep --job=cone20
    e4shared --run --job=cone20
    e4shared --post --job=cone20 --vtk-xml

If all of that worked successfully, you may view the final flow field using Paraview::

    paraview plot/cone20.pvd

The image below shows contours of pressure in the flow domain.

![Pressure contours for supersonic flow over a cone](/images/cone20-p-contour.png)


## Updating your installed version of Eilmer

Because Eilmer, L1d and the toolkit functions are being actively developed,
we will make frequent changes
and additions to the source code and the examples.
To update your copy of the repository,
move into any directory within the working tree and pull any new revisions into your local copy.

    cd dgd/src/eilmer
    make clean
    git pull -v

Then, to build a refreshed copy of Eilmer,
use `make` to coordinate the compiling and installing as before:

    make install


## Where to go from here

<!---
From here, you might like to look at more [Examples](/examples) or dive
into the [Documentation](/documentation).
-->
From here, you might like to dive
into the [Documentation]({{< relref documentation >}}).


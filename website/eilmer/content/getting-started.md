---
title: "Getting Started"
date: 2017-07-18T22:19:47+10:00
draft: false
menu: 
   main:
     weight: 2
---

## Prerequisites: background

This section is about you, the user of Eilmer. 
We assume that your mathematics, science or engineering background
adequately prepares you for computational fluid dynamics (CFD) analysis.
In particular, we assume that you have a working knowledge of geometry, calculus, mechanics,
and thermo-fluid-dynamics, at least to a second- or third-year university level.
With the Eilmer code, we try to make the analysis of compressible, reacting flow accessible
and reliable; we cannot make it trivial.

## Prerequisites: software

The core solver and its modules are mainly written in the
D programming language for speed and the benefits of compile-time checking. 
The pre- and post-processing modes make use of the Lua scripting language
so that we get flexibility and convenient customization.
There is also some Ruby a little Tcl/Tk used in the automated testing scripts.

To run simulations, you need an executable version of the Eilmer program.
<!--- RJG, 2018-05-23: Commented out while tarball unavailable in lead-up
                       to 4.0 release
You may get this executable file in the pre-built (tar-ball) package
linked to from the [front page](/) of this web site.
Alternatively,
--->
You may build an executable version of the program from the source code, as described below.

<!--- See comment above

If building from source is your choice,
continue reading the remainder of this section.
If you will be happy running the pre-built version,
go back to the [front page](/) and follow the links there,
but then come back here to
continue reading from "Setting up for first-time run" below.
--->

## Prerequisites: building from source

Our main development environment is Linux but the programs can be deployed on
Linux, flavours of Unix such as MacOS-X, and MS-Windows.
The main requirement is a D language compiler.
The source code of the Lua interpreter is included in the Eilmer source code repository.

To build Eilmer, you will require:

  + A C compiler (GNU compiler is a good option and standard on most systems)
  + A D compiler (The DMD compiler is a good choice, ldmd2 works equally well)
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
  + plotutils development package:
      + libplot-dev on Debian/Ubuntu/Mint
      + plotutils-devel on RedHat/CentOS/Fedora

Additionally, if you want to run the test suite, you will require:

  + Ruby package
  + TCL package
  + the Python sympy package

For viewing and plotting results, we recommend:

  + Paraview
  + Gnuplot

### Installing on Mac OSX
A guide to getting your development environment set up on 
Mac OSX is available [here](/mac-install).


## Getting the source code

The full source code for Eilmer4 and a set of examples can be found in a public repository
on [bitbucket](https://bitbucket.org/cfcfd/dgd-git).
To get your own copy, use the git revision control client to clone the repository 
with something like the following command:

    git clone https://bitbucket.org/cfcfd/dgd-git dgd

and within a couple of minutes, depending on the speed of your network connection,
you should have your own copy of the full source tree and the complete repository history.

## Installing

The default installation directory is `$HOME/dgdinst`.
To compile and install Eilmer, move into the eilmer source
area and use `make` to coordinate the compiling and installing:

    cd dgd/src/eilmer
    make install

If you are on a Mac, you'll need to give the `make` command an
extra hint:

    make PLATFORM=macosx install

## Setting up for first-time run

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


## Running your first simulation

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

Because Eilmer is being actively developed, we will make frequent changes
and additions to the source code and the examples.
To update your copy of the program and examples, move into the eilmer source
area and pull any new revisions into your local copy.
Then use `make` to coordinate the compiling and installing as before:

    cd dgd/src/eilmer
    make clean
    git pull -v
    make install


## Where to go from here

<!---
From here, you might like to look at more [Examples](/examples) or dive
into the [Documentation](/documentation).
-->
From here, you might like to dive
into the [Documentation](/documentation).


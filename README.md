# GDTk -- Gas Dynamics Toolkit
This repository hosts our collection of tools, 
written mainly in the D programming language, 
for gas dynamic simulations.
Our focus is on open source development to give a simple
access point for doing gas dynamics in research and teaching.


## Contents
* Toolkit Overview
* (Relatively) Quick Start
* Eilmer -- A compressible flow solver
    * Features
* Documentation
* License
* Contributors
* Chief Gardeners

### Toolkit Overview
The Gas Dynamics Toolkit is a collection of programs and functions for
computing the properties of high-speed gas flows.
Since the computational tools have been developed within the University of Queensland's
Centre for Hypersonics, there is clear bias toward predicting chemically-reacting gas flows,
as would be found in shock-tunnel and expansion-tube experiments.

The computational tools range from large-scale programs for the simulation
of gas flows in whole experimental facilites
and more limited functions for the evaluation of simple state-to-state transitions,
such as the jump in flow conditions across a simple shock wave.


### Who we are

The principal developers of these tools are Rowan Gollan, Kyle Damm, Nick Gibbons, and Peter Jacobs.
Substantial support has been given to the project by Anand Veeraragavan, Ingo Jahn, Vince Wheatley, and Fabian Zander.
There have been many more contributors to the project over the years,
including colleagues at other universities, students and visitors.

### Tools

The tools in the kit are:

* Eilmer,
  a simulation program for 2- and 3-dimensional gas flows.
* L1d,
  a program for the end-to-end simulation of free-piston-driven shock tunnels
  and expansion tubes.
* Pitot,
  a program using state-to-state calculations to estimate test flow conditions
  in impulse facilities.
* ESTCN,
  a state-to-state calculation program for estimating flow conditions in reflected-shock tunnels.
* NENZF1d,
  a program for estimating flow conditions in reflected-shock tunnels,
  when the test gas reaches temperatures high enough for chemical reactions to occur
  and when nonequilibrium chemistry effects are expected to be important.
* ceq, a calculator of thermochemical equilibrium gas compositions.
* Loadable Python libraries


### Prerequisite knowledge
This section is about you, the user of the Gas Dynamics Toolkit.
We assume that your mathematics, science or engineering background
adequately prepares you for computational fluid dynamics (CFD) analysis.
In particular, we assume that you have a working knowledge of geometry, calculus, mechanics,
and thermo-fluid-dynamics, at least to a second- or third-year university level.
With the Toolkit code, we try to make the analysis of compressible, reacting flow accessible
and reliable; we cannot make it trivial.

## Quick Start

This quick start guide introduces you to *Eilmer(,
and takes you through to running a simple simulation
of supersonic flow over a sharp cone.
It should take less than half an hour to do the full set up and run.

> **Note:** The section on installing prerequisites assumes a working knowledge
> on linux system configuration.
> If you find this section it too quick or too terse,
> you could try the gentler chapter on Installation.

### Install prerequisites
The main requirement is a D language compiler.
We recommend using the latest stable release of the LLVM D compiler.

To build Eilmer and other programs in the toolkit, you will require:

  + D compiler
      + Binary releases for the latest stable release of the LLVM D compiler (`ldc2` and `ldmd2`)
        may be found at: https://github.com/ldc-developers/ldc/releases .
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

### Getting the source code

The full source code for the toolkit programs, including a set of examples,
can be found in a public repository on [github](https://github.com/gdtk-uq/gdtk).
To get your own copy, use the git revision control client to clone the repository
with the following command:

    git clone https://github.com/gdtk-uq/gdtk.git gdtk 

and within a couple of minutes, depending on the speed of your network connection,
you should have your own copy of the full source tree and the complete repository history.

### Installing Eilmer

The default installation directory is `$HOME/gdtkinst`.
To compile and install *Eilmer*, move into the eilmer source
area and use `make` to coordinate the compiling and installing:

    cd gdtk/src/eilmer
    make install

If you are on a Mac, you'll need to give the `make` command an
extra hint:

    make PLATFORM=macosx install
    

### Configure environment

We'll assume you are happy using the default install area `$HOME/gdtkinst`.
The next step is to configure your environment to use Eilmer.
You will need to set the `DGD` variable to point to the top of the installation tree,
and the `DGD_REPO` variable to point to the top of the repository tree.
Note that the installation tree and repository tree are separate.
You then also need to set `$PATH`, `$DGD_LUA_PATH` and `$DGD_LUA_CPATH`
to point to the appropriate places.
Some example lines from a `.bashrc` file are:

    export DGD=$HOME/gdtkinst
    export DGD_REPO=$HOME/gdtk
    export PATH=$PATH:$DGD/bin
    export DGD_LUA_PATH=$DGD/lib/?.lua
    export DGD_LUA_CPATH=$DGD/lib/?.so
    
Remember to refresh your current shell (or log out and log back in) so
that your newly configured environment is available.

### Running your first *Eilmer* simulation

To test that everything has worked, you can exercise the flow
solver to simulate the supersonic flow over a 20-deg cone.

    cd ~
    cd gdtk/examples/eilmer/2D/sharp-cone-20-degrees/sg
    prep-gas ideal-air.inp ideal-air-gas-model.lua
    e4shared --prep --job=cone20
    e4shared --run --job=cone20
    e4shared --post --job=cone20 --vtk-xml

If all of that worked successfully, you may view the final flow field using Paraview::

    paraview plot/cone20.pvd


## Eilmer -- A compressible flow solver
Presently, the principal code in this collection is 
the *Eilmer* simulation code for 2D and 3D gas dynamic flows 
that may involve chemical reactions.
It is a research/education code and, 
with its built-in grid generation capabilities, 
is suitable for the exploration of flows 
where the bounding geometry is not too complex.

### Features:
* Eulerian/Lagrangian description of the flow (finite-volume, 2D axisymmetric or 3D).
* Transient, time-accurate, optionally implicit updates for steady flow.
* Shock capturing plus shock fitting boundary.
* Multiple block, structured and unstructured grids.
* Parallel computation in a shared-memory context.
* High-temperature nonequilibrium thermochemistry.
* GPU acceleration of the finite-rate chemistry.
* A selection of thermochemical models.
* Rotating frame of reference for turbomachine modelling.
* Turbulence models: Spalart-Allmaras and k-omega
* Conjugate heat transfer to solid surfaces and heat flow within 
  solid objects.
* Adjoint solver for design optimisation.
* MHD simulation for a single-fluid plasma (a work in progress).
* Import of GridPro structured grids and SU2 unstructured grids 
  for complex flow geometries.

We have structured the code as a programmable program, 
with a user-supplied input script (written in Lua) 
providing the configuration for any particular simulation exercise.
Our target audience is the *advanced* student of gas dynamics,
possibly an undergraduate student of engineering but, more likely,
a postgraduate student or academic colleague 
wanting to simulate gas flows as part of their study.

## Documentation
The documentation for users of the code is in a set of PDF reports
at the GDTk web site <http://gdtk.uqcloud.net>
Presently there are user guides for the main simulation code,
the geometry package and the gas model package.
More documents will appear as they are completed.

For those brave souls prepared to dive into the use and extension of the
code, there are examples provided as well as the source code itself.
For the short term, with the code in a period of rapid development,
we expect that users of the code will be mainly our
students and academic colleagues who can talk to us directly.

## License
For the source code, we use the GNU General Public License 3.
Please see the file `gpl.txt`.
For the documentation, we use the Creative Commons 
Attribution-ShareAlike 4.0 International License.

## Contributors
This code is the product of many people. 
Please see the file `contributors.txt` for some details of 
who has written what part of the code. 
The commit history is place to go to see further details. 

## Chief Gardeners
Peter Jacobs and Rowan Gollan, 2015-08-31 -- 2021-11-10


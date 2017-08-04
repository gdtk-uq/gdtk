# dgd -- D Gas Dynamics
This repository hosts our collection of tools, 
written mainly in the D programming language, 
to simulate compressible flows. 
Our focus is open source code that is designed to be a simple access 
point for teaching students about CFD codes.

## Contents
* Eilmer -- A compressible flow solver
    * Features
    * Documentation
    * Build Requirements
    * Running the program
* License
* Contributors
* Chief Gardeners

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
* Dense-gas thermodynamic models and rotating frames of reference for
  turbomachine modelling.
* Turbulence models.
* Conjugate heat transfer to solid surfaces and heat flow within 
  solid objects.
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

### Documentation
The documentation for users of the code is in a set of PDF reports
at the Eilmer web site http://cfcfd.mechmining.uq.edu.au/eilmer/
Presently there are user guides for the main simulation code,
the geometry package and the gas model package.
More documents will appear as they are completed.

For those brave souls prepared to dive into the use and extension of the
code, there are examples provided as well as the source code itself.
For the short term, with the code in a period of rapid development,
we expect that users of the code will be mainly our
students and academic colleagues who can talk to us directly.

### Build Requirements
Once you have cloned this repository, all that is required is 
a Linux environment with a fairly recent D compiler for building the
main code, a C compiler for building the Lua interpreter, and
a Fortran-2003 compiler for building some of the thermochemical models.

You may use the reference DMD compiler or the GDC or LDMD2 compilers.
We have been developing with the recent DMD32 and DMD64 compilers, 
from version 2.073.
For the C and Fortran compilers, we use gcc and gfortran.

Beyond having a standard Linux system with recent compilers,
the build of the Lua interpreter requires the development packages for
libreadline and libncurses.

Going into the `dgd/src/eilmer` directory you will find a single `makefile`
that allows the build to proceed with the command `make install`.
The executable program and supporting files will be installed into the 
directory `$HOME/dgdinst/` by default.

### Running the program
For running the program, environment variables 
may be set for the bash shell with:

    export DGD_REPO=${HOME}/dgd
    export DGD=${HOME}/dgdinst
    export PATH=${PATH}:${DGD}/bin
    export DGD_LUA_PATH=${DGD}/lib/?.lua
    export DGD_LUA_CPATH=${DGD}/lib/?.so

Setting DGD_REPO may be handy if you have cloned your copy of the
repository to somewhere other than `$HOME/dgd/`.

The actual running of a simulation is done in stages.

1. Prepare the grids and initial flow configuration by interpreting your
   Lua input script.
2. Run the main simulation program, starting from the prepared initial 
   flow state and allowing the flow field to develop in time, writing
   the resulting flow field states at particular times.
3. Postprocess the simulation data to extract particular data of interest.

Of course, this description is too superficial to actually expect that 
you will be able to run a simulation with no further instruction.
If you are keen, it's now time to look at the PDF documents 
mentioned above and then look at specific examples provided in the 
`dgd/examples/eilmer/` directory.

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
Peter Jacobs and Rowan Gollan, 2015-08-31 -- 2017-06-05


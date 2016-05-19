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
is suitable for the exploration of flows where the bounding geometry 
is not too complex.

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
* Coupling to radiation and ablation codes for aeroshell flows.
* Conjugate heat transfer to solid surfaces and heat flow within 
  solid objects.
* MHD simulation for a single-fluid plasma.
* Import of GridPro grids for complex flow geometries.

We have structured the code as a programmable program, 
with a user-supplied input script (written in Lua) 
providing the configuration for any particular simulation exercise.
Our target audience is the *advanced* student of gas dynamics,
possibly an undergraduate student of engineering but, more likely,
a postgraduate student or academic colleague 
wanting to simulate gas flows as part of their study.

### Documentation
Although there's not much documentation currently available for this
D-language code collection, the user guides from 
the [CFCFD](http://cfcfd.mechmining.uq.edu.au) project will
give some idea of the use of this new code.
In particular, 
the [Eilmer3 User Guide](http://cfcfd.mechmining.uq.edu.au/pdf/eilmer3-user-guide.pdf)
shows how to use the previous C++/Python code 
and [Theory Book](http://cfcfd.mechmining.uq.edu.au/pdf/eilmer3-theory-book.pdf),
describes how it was built.

For those brave souls prepared to dive into the use and extension of the
code, there are examples provided as well as the source code itself.
For the short term, with the code in a period of rapid development,
we expect that users of the code will be mainly our
students and academic colleagues who can talk to us directly.

### Build Requirements
Once you have cloned this repository, 
all that is required Linux environment with a fairly recent D compiler
and a C compiler (for building the Lua interpreter).
You may use the reference DMD compiler or the GDC or LDMD2 compilers.
We have been developing with the DMD32 and DMD64 compilers 
from versions 2.067 through 2.070.

Going into the `dgd/src/eilmer` directory you will find a single `makefile`
that allows the build to proceed with the command `make install`.
The executable program and supporting files will be installed into the 
directory `$HOME/dgdinst/` by default.

### Running the program
For running the program, environment variables 
may be set for the bash shell with:

    export DGD=$HOME/dgdinst
    export PATH=$PATH:$DGD/bin
    export DGD_LUA_PATH=$DGD/lib/?.lua
    export DGD_LUA_CPATH=$DGD/lib/?.so

The actual running of a simulation is done in stages.

1. Prepare the grids and initial flow configuration by interpreting your
   Lua input script.
2. Run the main simulation program, starting from the prepared initial 
   flow state and allowing the flow field to develop in time, writing
   the resulting flow field states at particular times.
3. Postprocess the simulation data to extract particular data of interest.

Of course, this description is too superficial to actually expect that 
you will be able to run a simulation with no further instruction.
If you are keen, it's now time to look at the Eilmer3 documents 
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
Peter Jacobs and Rowan Gollan, 2015-08-31 -- 2016-05-19


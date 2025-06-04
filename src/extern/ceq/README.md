# equilibrium-c -- A lightweight, modern library for equilibrium chemistry calculations

- Author: Nick Gibbons (n.gibbons(at)uq.edu.au)

- References:

    "equilibrium-c: A Lightweight Modern Equilibrium Chemistry Calculator for Hypersonic Flow Applications"\
    Preprint available on [arxiv.org](https://arxiv.org/abs/2412.07166), 2023\
    Nicholas N. Gibbons

    "Computer Program for Calculation of Complex Equilibrium Compositions and Applications"\
    NASA Reference Publication 1311, October 1995\
    Sanford Gordon and Bonnie J. McBride

    "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species"\
    NASA/TP - 2002-211556, September 2002\
    Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon

- Build Requirements

    + python3
    + python3-numpy
    + gcc

- Build Instructions (Linux)

    To build and install in '$HOME/eqcinst' type:\
    $ cd source\
    $ make all\
    $ make install


    Installation in a custom directory is accomplished by the INSTALL_DIR variable:\
    $ make install INSTALL_DIR='/path/to/install/dir'


    Once installed, add the location to your PYTHONPATH in your .bashrc file:\
    export PYTHONPATH=${PYTHONPATH}:${HOME}/eqcinst

- Build Instructions (Windows)

    Please click [here](./build_for_windows.md)


- Use Instructions

    The code is accessed through a python script that handles reading of data, memory management, and initialisation. All of this data is then passed to a set of c routines that do most of the actual work. See examples in the "examples" directory.

```python
from numpy import array
import eqc

species = ['CO2', 'CO', 'O2']
eq = eqc.EqCalculator(species)

X0 = array([1.0, 0.0, 0.0])
X1 = eq.pt(p=10e3, T=2500.0, Xs0=X0)
print(X1)

>>> array([ 0.65897178,  0.22735215,  0.11367607])
```

- Licensing

    This program is MIT licensed, which allows you do to almost anything with it except pretend that you wrote it or mess with the license itself. See the mit.txt file for details.


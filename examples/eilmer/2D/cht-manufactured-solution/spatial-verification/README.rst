Notes on running the Method of Manufactured Solutions coupled domain case
=========================================================================
:Author: Rowan J. Gollan
:Date: 2015-08-30

This test case is used to verify that coupled gas/solid solver is giving
an observed order of spatial accuracy of approximately 2.0. The layout
and use of this verification exercise is similar to the fluid-only
manufactured solutions case in ``examples/eilmer4/2D/manufactured-solution``.

How to run the verification case
--------------------------------
A Python script is used to coordinate the execution of this test case
on a number of grids in succession. The chosen grid sizes are specified
in a configuration file. There is an example configuration file in this
directory that will run the following grids: 8x12, 16x24, 32x48 and 64x96.
When specifying the grid sizes, the cell count for the ``i`` direction
is given. The cell count in the ``j`` direction is one and half times the
number of cells in the ``i`` direction since the domain is one and half
times as high as it is wide. This gives square cells. 
To run this example, do::

  > ./run-verification-test.py example-config.py

During execution, the script will run each of the grids in a separate
subdirectory. Those names follow the format of ``nicellsxnjcells``.

At the end of a succesful run, a set of files are produced that give the
L2 and Linf norms based on temperature error and the observed order of
spatial accuracy. These will be sitting in the directory where the script
was started.

Configuring the cases to run
----------------------------
The ``run-verification-script.py`` requires one argument and that is the
name of a configuration file. The configuration file is a Python file that
contains the specifics of the case to run. An example file is shown below
and details on each of the parameters are listed afterwards.

::

  ncellsList = [8, 16, 32, 64]

ncellsList : list
  This list gives the sequence of grids. Each grid is designated
  by a single integer. A value of 8 indicates that a grid with
  8 x 12 cells will be used. Note the value corresponds to the number
  of cells in the ``i``direction. The list of values should be given
  in increasing order. It is typical, but not mandatory, to use
  use a refinement ratio of 2 by giving a doubling of the number
  of cells for each successive grid.




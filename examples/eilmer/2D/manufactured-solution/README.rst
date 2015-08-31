Notes on running the Method of Manufactured Solutions case
==========================================================
:Author: Rowan J. Gollan and Peter A. Jacobs
:Date: 2015-08-07

This test case is taken from the paper by Roy et al. (2004).
We used this test case with success as a verification exercise
for Eilmer3, reported in Gollan and Jacobs (2013). In this
directory, you will find updated scripts to run this case
with Eilmer4.

How to run the verification case
--------------------------------
Since this is a verification case, we are interested in extracting
an observed order of accuracy. The observed order of accuracy test
is widely regarded as the most rigorous form of verification. In
this test, we check that the spatial order of accuracy of the code
by comparison to an exact (manufactured) solution. This "observed"
order of accuracy is then compared against the formal order of
accuracy of our chosen numerical methods. If the values for order
of accuracy agree (in the limit of diminishing cell size), then
that is a strong indication that the implemented numerical solution
technique is correct. As such, we need to run this case on several
grids in  succession so that we can extract the observed order of
accuracy. To help with that task, there is a script that can be run.
The script coordinates the running of several grids. It determines
what cases you want run based on a user-supplied configuration file.
There is more detail on the configuration file in the next section.
For now, you could try the example configuration file that runs the
Euler case on a single processor by executing the script at the
command line:

  > ./run-verification-test.py example-config.py

During execution, the script will run each of the grids in a
separate subdirectory. The subdirectory names were graciously 
suppliedby the Department of the Bleeding Obvious. Those names
follow the format of ``ncellsxncells``.

At the end of a succesful run, a set of files are produced that
give the L2 and Linf norms based on density error and the observed
order of spatial accuracy. These will be sitting in the directory
where the script was started.

Configuring the cases to run
----------------------------
The ``run-verification-script.py`` requires one argument and that is
the name of a configuration file. The configuration file is a Python
file that contains the specifics of the case to run. An example file
is shown below and details on each of the parameters are listed
afterwards.

::

  case = 1
  ncellsList = [8, 16, 32, 64]
  fluxCalc = 'ausmdv'
  blocking = 'single'
  threading = 'multi'


case : integer
  The case is used to select between the Euler and Navier-Stokes
  versions of the manufactured solution, and whether or not
  a scaling function is applied. The scaling function may be used
  to de-emphasise any corner and edge effects and concentrate
  the error on the interior of the domain. Allowable values are:
  1 -- Euler case; 2 -- Navier-Stokes case; 3 -- Euler case with scaling;
  and 4 -- Navier-Stokes case with scaling.

ncellsList : list
  This list gives the sequence of grids. Each grid is designated
  by a single integer. A value of 8 indicates that a grid with
  8 x 8 cells will be used. The list of values should be given
  in increasing order. It is typical, but not mandatory, to use
  use a refinement ratio of 2 by giving a doubling of the number
  of cells for each successive grid.

fluxCalc : string
  This value is used to select which flux calculator is tested.
  The available flux calculators are (at the time of writing):
  'ausmdv', 'efm', 'ausm_plus_up', 'adaptive', and 'hlle'.

blocking : 'string'
  This value indicates how the domain will be split with the allowable
  strings being 'single' or 'multi'. The 'single' selection denotes that
  a single block is used to cover the domain. This is useful to check the
  accuracy absent any compilcations from internal corner connections that
  occur in multiblock simulations. The 'multi' option will split the
  domain into four blocks, two in the ``i`` direction and two in the ``j``
  direction.

threading: 'string'
  This value indicates whether the simulation will be run on a single
  thread or on multiple threads. A value of 'single' will restrict Eilmer4
  to using a single thread by setting ``--max-cpus=1`` when executing the
  simulation program. When 'multi' is selected, Eilmer4 will use its
  shared memory, multiple thread mode for execution. The number of threads
  used is the smaller of the number of blocks or the total number of
  execution threads (think cores, on non-hyperthreaded machines) available
  on the machine. That is: ``min(no. blocks, max no. execution threads)``.
  Note that running the Navier-Stokes case on a single thread gets time
  consuming very quickly as the number of cells are increased.

  The threading option is mostly used by the developers to check that the
  code gives identical results when run serially or in parallel. If things
  are working correctly, there should be no difference in observed order of
  spatial accuracy due to running in either serial or parallel mode.

References
----------

::
  Gollan, R.J. and Jacobs, P.A. (2013)
  About the formulation, verification and validation of the hypersonic flow solver Eilmer.
  International Journal for Numerical Methods in Fluids, 73:19--57, DOI: 10.1002/fld.3790.

  Roy, C.J., Nelson C.C., Smith T.M., and Ober C.C. (2004)
  Verification of Euler/Navier-Stokes codes using the method of manufactured solutions.
  International Journal for Numerical Methods in Fluids, 44:599–620.




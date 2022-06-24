Notes on the parallel-scaling example
=============================================
:Author: Nick Gibbons (n.gibbons@uq.edu.au)
:Date: 2022-06-24

This example is a bit unusual. The main directory, full-test-case, is a
simulation of a 15 degree blunt cone with a 5mm nose radius, in Mach 5 air at
300K and 10 kPa. The mesh is parametrically generated in full 3D, and can be
scaled up or down in fidelity to add as many cells as the user desires.

The actual purpose of this example can be found in the other directories, e4mpi
and nk, which contain modified copies of the sphere cone calculation that can
be used to benchmark the code's parallel performace. This is accomplished by
running multiple calculations of the same flow with different numbers of
processors working on each one. The iterations per second of the simulations
with more processors should be faster than those with fewer, but with some
degradation as the number is increased. Running the example generates some
empirical data to estimate the magnitude of this degradation.

To perform your own scaling tests, simply alter the numbers in parameters.txt,
which sets the schedule of the test simulations that will be run. For example,
the default value is number_of_processors=1,2,4,8,12,16; which runs six
simulations in total, one with 1 core, one with 2, one with 4 etc. Using run.sh
will then execute a series of python scripts that create the tests, run them,
and process the results.

To run with a larger number of cores, perhaps for a supercomputer, you may want
to try editing the gengrid.lua file to add more cells to the problem. You
should be able to do this using the N_REFINE parameter, although you may also
have to change a_factor on line 165 to keep the cluster functions happy.

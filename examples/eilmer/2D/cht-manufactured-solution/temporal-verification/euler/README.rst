Notes on running the Method of Manufactured Solutions coupled domain case
=========================================================================
:Author: Anand V. and Rowan G.
:Date: 2015-08-30

This test case is used to verify that the coupled gas/solid solver is giving
the expected order of temporal accuracy. This particular case tests
the Euler method, in which we expect an order of temporal accuracty of 1.

How to run the verification case
--------------------------------
A Python script is used to coordinate the execution of this test case
for a number of timestep sizes. The chosen timestep sizes are specified
in a configuration file. There is an example configuration file in this
directory that will run the following timesteps: 1.0e-6, 0.5e-6, 0.25e-6, 0.125e-6.
To run this example, do::

  > ./run-verification-test.py example-config.py

During execution, the script will run each of the timestep cases in a separate
subdirectory. Those names follow the format of ``dt-<CASE-ID>``.

At the end of a succesful run, a set of files are produced that give the
L2 and Linf norms based on temperature error and the observed order of
temporal accuracy. These will be sitting in the directory where the script
was started.

Configuring the cases to run
----------------------------
The ``run-verification-script.py`` requires one argument and that is the
name of a configuration file. The configuration file is a Python file that
contains the specifics of the case to run. An example file is shown below
and details on each of the parameters are listed afterwards.

::

  dtList = [1.0e-6, 0.5e-6, 0.25e-6, 0.125e-6]

dtList : list
  This list gives the sequence of timesteps to use. Each timestep case is designated
  by a float value. It represents a timestep value in seconds.
  It is typical, but not mandatory, to halve the timestep for each successive case.




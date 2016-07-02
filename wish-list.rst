Wish list for DGD code collection
=================================
:Authors: Rowan J. Gollan and Peter A. Jacobs
:Date: 2016-02-05

This document contains a wish list of desired items and code
to support the DGD project. This could be useful starting
point for new developers who wish to contribute to the project.
You're encouraged to pick a task and work on it. It would
be best to correspond with one of the lead developers first
so that we can send you down the right garden path, and to
be sure that the item is still part of the wish list.

The rest of this document is basically a large TODO list.
It has been broken into sections based on theme. Not all tasks
are strictly code development. Some are related to supporting
in the project in other ways such a webpage curation, documentation
and illustration generation.


Marketing
---------

Webpage development
  Along with our entry on bitbucket, we should have a set of webpages
  hosted on the UQ servers that advertise our project. I think that
  using a Sphinx would be a good way to manage this, but I'm not tied
  to that tool specifically.

Eilmer
------

Source code development
^^^^^^^^^^^^^^^^^^^^^^^
+ Update spatial derivative evaluation in solid domain to emulate that
  in the gas domain.
+ Extend gas/solid block connections to all orientations in 2D.
+ Implement and test a moving wall boundary condition.
+ Implement and test k-omega model.
+ Implement and test two-temperature modelling.
+ Implement cell search in 3D for ``setHistoryPoint``.
+ Implement diffusion models, of particular interest are:
  - Fickian diffusion with constant Lewis numbers.
  - Fickian diffusion with mixture-averaged diffusion coefficients.
  - Armaly-Sutton diffusion model.

User interface
^^^^^^^^^^^^^^
+ A rendering tool/approach to build a static image of the computational
  domain in 2D and 3D.
+ A rendering tool to build a dynamic view of the computational domain.


Verification
^^^^^^^^^^^^
+ Port over oblique detonation wave example from eilmer3.
+ Use the Method of Manufactured Solutions to test 3D implementation.
+ Investigate means to verify boundary conditions.

Documentation
^^^^^^^^^^^^^
+ Build a video tutorial.
+ Build and maintain an FAQ.
+ Proof reading of user guide.
+ Proof reading of examples.
+ Addition of new examples.

Gas package
-----------

Source code
^^^^^^^^^^^
+ Examine Lua-wrappers with a view to improving the consistency of the interaction
  with users.
+ Increase unit-test coverage.

Data
^^^^
+ Build a script to drive the ``species-data-converter`` to bring over all
  the species required for methane combustion.

Documentation
^^^^^^^^^^^^^
+ Proof-read *gas* package user guide.
+ Build an example of using ``gas-calc`` to do a thermo cycle calculation.

Kinetics package
----------------

Source code
^^^^^^^^^^^
+ Work on the Lua-wrappers to expose the kinetics routines.

Reaction schemes
^^^^^^^^^^^^^^^^
+ Port over air reaction schemes.
+ Port over hydrogen combustion schemes.
+ Port over methane combustion schemes.

Documentation
^^^^^^^^^^^^^
+ Build examples of using the *kinetics* pacakge as a stand-alone calculator.





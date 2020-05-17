---
weight: 10
title: Gas dynamic tools
---

# The gas-dynamic tools
The principal tools in the kit are:
+ [Eilmer]({{< relref "eilmer/" >}}),
  a simulation program for 2- and 3-dimensional gas flows.
+ [L1d]({{< relref "l1d/" >}}),
  a program for the end-to-end simulation of free-piston-driven shock tunnels
  and expansion tubes.
+ [foamMesh]({{< relref "foammesh" >}}),
  A nice parametric grid-generation tool for OpenFOAM simulations.
+ [IMOC]({{< relref "moc.md" >}}),
  an interactive program for the isentropic method for characteristics
  for 2-dimensional ideal gas flows.
+ A library of functions for the calculation of simple gas dynamic processes
  such as shocks and expansions in the presence of chemical reactions.
  This library may be used from within your Lua scripts,
  or a version of it may loaded into
  a [Python]({{< relref "python-library.md" >}}) or
  Ruby interpreter.
  The Poshax program is an example of a Python program that makes use of the functions
  from the library.



---
title: "GDTk: The Gas Dynamic Toolkit"
---

# GDTk: The Gas Dynamics Toolkit

The Gas Dynamics Toolkit is a collection of programs and functions for
computing the properties of high-speed gas flows.
Since the computational tools have been developed within the University of Queensalnd's
Centre for Hypersonics, there is clear bias toward predicting chemically-reacting gas flows,
as would be found in shock-tunnel and expansion-tube experiments.

The computational tools range from large-scale programs for the simulation
of gas flows in whole experimental facilites
and more limited functions for the evaluation of simple state-to-state transitions,
such as the jump in flow conditions across a simple shock wave.

The principal developers of these tools are Peter Jacobs, Rowan Gollan,
Ingo Jahn, Anand Veeraragavan and Vince Wheatley.
There have been many more contributors to the project over the years,
including colleagues at other universities, students and visitors.
A list of some of those contributors is available [here]({{< relref contributors >}}).

The principal tools in the kit are:
+ [Eilmer]({{< relref "docs/tools/eilmer/" >}}),
  a simulation program for 2- and 3-dimensional gas flows.
+ [L1d]({{< relref "docs/tools/l1d/" >}}),
  a program for the end-to-end simulation of free-piston-driven shock tunnels
  and expansion tubes.
+ [IMOC]({{< relref "docs/tools/python-library/moc.md" >}}),
  an interactive program for the isentropic method for characteristics
  for 2-dimensional ideal gas flows.
+ A library of functions for the calculation of simple gas dynamic processes
  such as shocks and expansions in the presence of chemical reactions.
  This library may be used from within your Lua scripts,
  or a version of it may loaded into
  a [Python]({{< relref "docs/tools/python-library/" >}}) or
  Ruby interpreter.
  The Poshax program is an example of a Python program that makes use of the functions
  from the library.


## How to get started
The complete repository is available on bitbucket: https://bitbucket.org/cfcfd/dgd-git

It is worth browsing the [Getting started]({{< relref getting-started >}}) introduction
and trying the starting exercise described there.


## License
The codes and documentation in the this toolkit are freely available.
+ For the source code within the toolkit, we use the GNU General Public License 3.
Please see the file `gpl.txt` in the source tree.
+ For the documentation, such as this website and the user guides,
we use the Creative Commons Attribution-ShareAlike 4.0 International License.


## How to cite
We hope that, by using Eilmer and the other tools,
you are able to produce some high quality simulations that aid your work.
When it comes time to report the results of your simulations to others,
we ask that you acknowledge our work by citing our papers on the codes.


## How can you contribute
+ By using the code and giving us feedback on your simulations
  and results.
+ By giving us feedback on the documentation.
+ By contributing source code.



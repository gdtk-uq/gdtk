---
title: "About the toolkit"
description: "GDTk"
lead: "GDTk is built by the compressible flow CFD group at The University of Queensland"
date: 2020-10-06T08:48:57+00:00
lastmod: 2020-10-06T08:48:57+00:00
draft: false
images: []
menu:
  docs:
    parent: "introduction"
weight: 10
toc: true
---

The Gas Dynamics Toolkit is a collection of programs and functions for
computing the properties of high-speed gas flows.
Since the computational tools have been developed within the University of Queensland's
Centre for Hypersonics, there is clear bias toward predicting chemically-reacting gas flows,
as would be found in shock-tunnel and expansion-tube experiments.

The computational tools range from large-scale programs for the simulation
of gas flows in whole experimental facilites
and more limited functions for the evaluation of simple state-to-state transitions,
such as the jump in flow conditions across a simple shock wave.


## Who we are

The [principal developers]({{< relref "dev-team-and-contributors" >}}) of
these tools are Rowan Gollan, Kyle Damm, Nick Gibbons, Daryl Bond and Peter Jacobs.
Substantial support has been given to the project by Anand Veeraragavan, Ingo Jahn, and Vince Wheatley.
There have been many more contributors to the project over the years,
including colleagues at other universities, students and visitors.
A list of those contributors is available [here.]({{< relref "dev-team-and-contributors#contributors" >}})

## Tools

The tools in the kit are:
+ [Eilmer]({{< relref "docs/eilmer/about" >}}),
  a simulation program for 2- and 3-dimensional gas flows.
+ [L1d]({{< relref "docs/l1d/about" >}}),
  a program for the end-to-end simulation of free-piston-driven shock tunnels
  and expansion tubes.
+ [Pitot]({{< relref "docs/pitot/about" >}}),
  a program using state-to-state calculations to estimate test flow conditions
  in impulse facilities.
+ [ESTCN]({{< relref "docs/estcn/estcn-manual-for-hugo" >}}),
  a state-to-state calculation program for estimating flow conditions in reflected-shock tunnels.
+ [NENZF1d]({{< relref "docs/nenzf1d/nenzf1d-manual-for-hugo" >}}),
  a program for estimating flow conditions in reflected-shock tunnels,
  when the test gas reaches temperatures high enough for chemical reactions to occur
  and when nonequilibrium chemistry effects are expected to be important.
+ [Loadable Python libraries]({{< relref "docs/python/loadable-library" >}})


## Prerequisite knowledge
This section is about you, the user of the Gas Dynamics Toolkit.
We assume that your mathematics, science or engineering background
adequately prepares you for computational fluid dynamics (CFD) analysis.
In particular, we assume that you have a working knowledge of geometry, calculus, mechanics,
and thermo-fluid-dynamics, at least to a second- or third-year university level.
With the Toolkit code, we try to make the analysis of compressible, reacting flow accessible
and reliable; we cannot make it trivial.

## How to get started
The complete repository is available on githubt: <https://github.com/gdtk-uq/gdtk>

If you are reasonably proficient with Linux (or Mac) at the command line,
then try the [Quick Start guide]({{< relref "quick-start" >}}).

Otherwise check out the [Prerequisite Software]({{< relref "docs/getting-started/prerequisites" >}})
and [Linux Install Guide]({{< relref "docs/getting-started/installation-locally" >}}) for more detailed
instructions.

## License
The codes and documentation in the this toolkit are freely available.
+ For the source code within the toolkit, we use the GNU General Public License 3.
Please see the file `gpl.txt` in the source tree.
+ For the documentation, such as this website and the user guides,
we use the Creative Commons Attribution-ShareAlike 4.0 International License.


## How can you contribute
+ By using the code and giving us feedback on your simulations
  and results.
+ By giving us feedback on the documentation.
+ By contributing source code.

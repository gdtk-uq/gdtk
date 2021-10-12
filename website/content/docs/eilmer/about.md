---
title: "Eilmer"
description: "About Eilmer"
lead: ""
date: 2020-10-06T08:48:57+00:00
lastmod: 2020-10-06T08:48:57+00:00
draft: false
images: []
menu:
  docs:
    parent: "eilmer"
weight: 10
toc: true
---

<span style="font-size:larger;">
<i>(pronounced: /&#603;lm&#601;/)</i>
</span>

{{< figure src="/images/ramp-in-3D-by-kad.png" width="100%" alt="scramjet flowpath" caption="<em>Simulation of flow over a 24° compression corner as in the experiments performed by Holden (1970). The simulation is performed in 3D as this captures the spillage of gas over side of the ramp and results in a better comparison to experimental results for heat transfer.</em>" >}}

## What does Eilmer do?

The Eilmer code is a program for the numerical simulation of transient,
compressible gas flows in two and three dimensions.
This program answers the "What if ... ?" type of question
where you will set up a flow situation by
defining the spatial domain in which the gas moves,
set an initial gas state throughout this domain,
specify boundary condition constraints to the edges of the domain and
then let the gas flow evolve according to the rules of gas dynamics.

Eilmer began life as a tool to simulate and aid in the of design shock tunnels and expansion tubes.
It has also been applied to the analysis of the experiments undertaken
in shock tunnels and expansion tubes and has being used for the simulation and design of
hypersonic inlets, turbomachinery and microcombustors.

## Features

+ 2D/3D compressible flow simulation
+ Gas models include ideal, thermally perfect, equilibrium, multi-temperature.
+ Finite-rate chemistry
+ Inviscid, laminar, turbulent (k-omega) flows.
+ Solid domains with conjugate heat transfer in 2D.
+ User-controlled moving grid capability.
+ Shock-fitting method for 2D geometries.
+ A rotating frame of reference for turbomachine modelling.
+ Transient, time-accurate, using explicit Euler, predictor-corrector and Runge-Kutta updates.
+ Steady-state solver using the Newton-Krylov approach.
+ Parallel computation using shared memory or MPI.
+ Multiple block, structured and unstructured grids.
+ Native grid generation and import capability.
+ Unstructured-mesh partitioning via Metis.

## Documentation

You should probably read the awesome documentation that goes with Eilmer.


---
weight: 10
---

# Eilmer: compressible flow simulation program

![](/images/kiock-mach-field.png)
Transonic flow through a plane turbine cascade (Kiock et al., 1986).
Simulation set up by Peter Blyton for Eilmer3 in 2011.
This calculation with Eilmer4 in 2020.
Visualization with Paraview.

## What does Eilmer do

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

- [Flow-Solver User Guide (PDF)](/pdfs/eilmer-user-guide.pdf)
- [Geometry Package User Guide (PDF)](/pdfs/geometry-user-guide.pdf)
- [Gas Package User Guide (PDF)](/pdfs/gas-user-guide.pdf)
- [Reacting Gas Thermochemistry with the thermally-perfect-gas model (PDF)](/pdfs/reacting-gas-guide.pdf)
- [Reference Manual (HTML)](/html/eilmer-reference-manual.html)
- [Catalog of examples](/html/eilmer-catalog-of-examples.html)



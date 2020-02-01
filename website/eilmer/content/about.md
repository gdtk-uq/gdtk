---
title: "About"
date: 2017-07-18T22:32:35+10:00
draft: false
menu: 
   main:
      weight: 1
---

## What does Eilmer do

The Eilmer code is a program for the numerical simulation of transient,
compressible gas flows in two and three dimensions.
This program answers the "What if ... ?" type of question 
where you will set up a flow situation by 
defining the spatial domain in which the gas moves,
set an initial gas state throughout this domain, 
specify boundary condition constraints to the edges of the domain and
then let the gas flow evolve according to the rules of gas dynamics.

## Who built Eilmer, and why

The principal developers of Eilmer are Peter Jacobs and Rowan Gollan,
both at the University of Queensland.
There have been many contributors to the project over the years.
The complete list of contributors is available [here](/contributors).

Eilmer began life as a tool to simulate and aid in the of design shock tunnels
and expansion tubes.
It has also been applied to the analysis of the experiments undertaken
in shock tunnels and expansion tubes.

More recently, Eilmer is being used in the simulation and design of
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

## How can you contribute

+ By using the code and giving us feedback on your simulations
  and results.
+ By giving us feedback on the documentation.
+ By contributing source code. The complete repository is available
  on bitbucket: https://bitbucket.org/cfcfd/dgd-git

## How to cite

We hope that by using Eilmer you are able to produce some high quality simulations
that aid your work.
When it comes time to report the results of your Eilmer simulations to others,
we ask that you acknowledge our work by citing our papers on the Eilmer code:

<p>
Jacobs P.A. and Gollan R.J. (2016)<br>
<cite>Implementation of a Compressible-Flow Simulation Code in the D Programming Language.</cite><br>
Advances of Computational Mechanics in Australia; 846:54-60<br>
 (DOI: 10.4028/www.scientific.net/AMM.846.54)</p>

<p>
Gollan R.J. and Jacobs P.A. (2013)<br>
<cite>About the formulation, verification and validation of the hypersonic flow solver Eilmer.</cite><br>
International Journal for Numerical Methods in Fluids; 73:19-57<br>
 (DOI: 10.1002/fld.3790)</p>




## License

For the source code, we use the GNU General Public License 3.
Please see the file `gpl.txt` in the source tree.

For the documentation, such as this website, we use the Creative Commons 
Attribution-ShareAlike 4.0 International License.

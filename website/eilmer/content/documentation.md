---
title: "Documentation"
date: 2017-07-18T22:35:06+10:00
draft: false
menu:
   main:
      weight: 15
---

The supporting documentation for **Eilmer** is available for download:

- [Eilmer4 User Guide](/pdfs/eilmer-user-guide.pdf)
- [Geometry Package User Guide](/pdfs/geometry-user-guide.pdf)
- [Gas Package User Guide](/pdfs/gas-user-guide.pdf)

The user guide for **foamMesh** is also available:

- [foamMesh User Guide](/pdfs/foammesh-user-guide.pdf)

Reference Guides
----------------
There are Reference Manuals for boths users and developers:

- [Eilmer Reference Manual for Users, v4.0](/html/users-reference-manual.html)



Presentations and seminars
--------------------------

Over the years, we have presented various aspects on the formulation,
physical models, applications and use of Eilmer. Presentation slides
are available for download. Some slides contain user input tips, but please note
that the finer details of user input may change over time.
These slides are not updated, but rather serve as historical record of activities.
Instead, check the User Guide for the definitive user input.

- [Estimating parallel compute performance for Eilmer simulations](/pdfs/cfh-seminar-oct-2019.pdf), 03 October 2019.  
  RJG seminar to the Centre for Hypersonics research group on estimating
  parallel performance in the context of Eilmer simulations.
  The presentation also includes discussion of the available
  load-balancing techniques and examples of their use.
  A simple method of timing tests is presented that can be used
  to give a quick assessment of parallel performance.

- [Progress of the Eilmer 4 transient flow solvers](/pdfs/eilmer-talk-pj-aug-2019.pdf), 29 August 2019.  
  PAJ seminar to the Centre for Hypersonics research group, with some discussion on MPI, complex numbers,
  one-sided flux calculator and moving grids,
  complex 2D grids and an oscillating-jet injector simulation.
  
- [User-defined Run-time Customisation for Eilmer Simulations](/pdfs/cfh-seminar-jul-2018.pdf), 19 July 2018.    
  RJG seminar on user-defined customisation in Eilmer to the Centre for Hypersonics at The University of Queensland.

- [Eilmer 4.0: A toolkit for doing gas-dynamic calculations](/pdfs/eilmer-talk-pj-july-2018.pdf), 5 July 2018.  
  PAJ seminar on applications of Eilmer for gas dynamics to the Centre for Hypersonics at The University of Queensland.
  Includes some discussion on running MPI jobs, a state-to-state Brayton cycle calculation,
  the rebuild of Poshax as a Lua script,
  and simulation of the X3R with the step in its shock tube.

- [Computational Hypersonics Research at The University of Queensland](/pdfs/rjg-seminar-strathclyde-2018-06-08.pdf), 8 June 2018.  
  RJG seminar presented to Department of Mechanical & Aerospace Engineering at the University of Strathclyde, Glasgow

- [Compressible flow short course at National University of Singapore, Day 1](/pdfs/nus-short-course-on-eilmer-day-1.pdf), 9 May 2018  
  [Compressible flow short course at National University of Singapore, Day 2](/pdfs/nus-short-course-on-eilmer-day-2.pdf), 10 May 2018  
  RJG presented a 2-day short course on compressible flow CFD using Eilmer at the National University of Singapore.

- [Eilmer4: A tool for analysing hypersonic flows](/pdfs/eilmer4-talk-pj-april-2017.pdf), Prepared April 2017, presented June 2017.
  PAJ seminar on the story of Eilmer, going back to the ICASE days.
  On the technical side, we look at the GasModel and ThermochemicalReactor classes in Eilmer4.
  Example based on the Powers and Aslam oblique detonation wave problem.

- [Eilmer4: the next step in the UQ simulation codes for high-enthalpy flows](/pdfs/eilmer4-talk-pj-dec-2016.pdf), December 2016.
  Presentation to the Osney Hypersonics Symposium of the Eilmer4 code development.
  A little bit on shared-memory parallelism.
  Example of the sphere heat-transfer again.
  
- [Eilmer4: A Computational Tool for High-Enthalpy Flow Simulation](/pdfs/eilmer4-talk-pj-sep-2016.pdf), September 2016.
  Presentation to the Franco-Australian Symposium in Paris.
  Example of stagnation-point heating of a sphere, with the shock-fitting procedure in action.

- [Implementation of a compressible-flow simulation code in the D programming language](/pdfs/eilmer4-talk-pj-may-2016.pdf), May 2016.
  PAJ seminar on the development of the Eilmer4 code over the past year.
  Some history (for motivation), formulation and code structure to take advantage of easy shared-memory parallelism in D.
  Example of sharp-nosed body to show input script format.
  Parallel performance of the forward-facing step revisited.
  Blunt-nosed cone for a real hypersonics application.

- [Implementation of a compressible-flow simulation code in the D programming language](eilmer4-talk-nov-2015.pdf), November 2015.
  First real showing of the Eilmer4 code outside of the Centre for Hypersonics.

- [Eilmer4: What’s happening?](/pdfs/eilmer4-talk-may-2015.pdf), May 2015.
  PAJ seminar to Centre for Hypersonics research group on development of Eilmer4, one year in.
  An overview of the finite-volume formulation and code structure.
  Some history of the Eilmer series of codes to motivate the recent development in D.
  Example of sharp-nosed body to show input script format.
  Parallel performance of the forward-facing step.

- [Eilmer3: a CFD code and its validation (A quarter of a century spent “doing sums”.)](/pdfs/eilmer3-talk-june-2014.pdf), June 2014.
  PAJ seminar to Centre for Hypersonics research group on the story of the Eilmer simulation code.
  A potted history of our CFD tools for compressible flow.
  Some stats on the code itself.
  Verification and validation exercises: sharp-nosed projectile, sphere in non-reacting flow,
  Imperial College convex ramp and CUBRC cylinder with flare.
  Code features for turbomachinery calculations.
  Now, let's build Eilmer4 in D.

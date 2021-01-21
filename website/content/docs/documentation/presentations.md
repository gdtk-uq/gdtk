---
weight: 30
---

# Presentations and Seminars

Over the years, we have presented various aspects on the formulation,
physical models, applications and use of Eilmer. Presentation slides
are available for download. Some slides contain user input tips, but please note
that the finer details of user input may change over time.
These slides are not updated, but rather serve as historical record of activities.
Instead, check the User Guide for the definitive user input.

- [Towards an Eilmer Four Point Zero Release](/pdfs/cfh-seminar-jan-2021.pdf), 21 January 2021.  
  RJG seminar to the Centre for Hypersonics research group about the plans for a fixed release
  of Eilmer in early 2021.
  The talk captures the motivation behind a fixed release model, what Eilmer 4.0 will include,
  and discusses progress to date towards that release.

- [Calculation of Test Flow Conditions for Reflected-Shock Tunnels](/pdfs/eilmer-talk-pj-2020-dec.pdf), 03 December 2020.  
  PAJ seminar to the Centre for Hypersonics research group on a couple of ways to estimate
  test flow conditions in reflected shock tunnels.
  The first uses Eilmer to do an axisymmetric simulation of the nozzle expansion and
  the second is a simpler quasi-one-dimensional calculation.
  Both make use of the finite-rate chemistry module.
  
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

- [The Case for Open Source in Research _and what we're doing in Eilmer_](/html/case-for-open-source.html), 10 March 2016.  
  RJG seminar in the Centre for Hypersonics seminar series.
  This talk makes a case for the use of open source codes in research and then links that to the efforts
  to support reproducible research when using Eilmer.

- [Implementation of a compressible-flow simulation code in the D programming language](/pdfs/eilmer4-talk-nov-2015.pdf), November 2015.  
  The 2nd Australasian Conference on Computational Mechanics, held at QUT,
  was the first real showing of the Eilmer4 code outside of the Centre for Hypersonics.
  A bit of history of the codes, finite-volume formulation, structure of the new code and a couple of examples.
  The Zucrow and Hoffman sharp-nosed projectile example is shown in reasonable detail
  so the the sections of the scripted input file can be described.

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

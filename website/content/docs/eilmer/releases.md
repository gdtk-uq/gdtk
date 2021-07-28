---
title: "Releases"
description: "Eilmer Releases"
lead: ""
date: 2021-07-28
lastmod: 2020-07-28
draft: false
images: []
menu:
  docs:
    parent: "eilmer"
weight: 15
toc: true
---

The Eilmer project uses a rolling release model.
You can keep up with the latest features and bug fixes by tracking the development
in the bitbucket repository.
We recommend this for advanced users or those that require developmental and
experimental feature sets.

We also offer stable releases.
Our model of a stable release is a social contract rather than a software engineering contract.
That social contract is basically a list of the features we consider stable and well-tested for general use.
More importantly, it is the list of features for which documentation is available
and we are willing and ready to support with bug fixes.

Why do we have this model?
Eilmer is a simulation code under constant development.
Some features are bleeding-edge, some features are experimental, some features lack documentation,
and some features are still buggy.
If you want to use these experimental features, we request patience
and understanding on the part of the users.
It won't be possible for us to have complete, bug-free, documented code on
all interesting features --- but we don't necessarily want to hold you back from using them.
If you fall into this category, be sure to keep up with the repository code on a regular basis.

If you are new to the code or prefer to stick with a set of stable features,
then this leads us to the concept of a "release version".
For each release, there will be a set of features that we consider stable, well tested
and well documented.
We are confident about supporting users with these features and will try to jump on any
identified bugs as quickly as possible.

This page will document the set of features supported in each release.

## Version 4.0

To checkout this version:

    TO COMPLETE WHEN TAGGED

  

Capabilities/features supported in v4.0.

+ transient time-stepping
  + Euler
  + predictor-corrector
  + RK-3 variants
+ local time-stepping
+ grid capabilities
  + structured grids
  + unstructured grids
  + moving grid (user-defined motion and shock-fitting)
  + import GridPro format
  + import SU2 format
+ parallel execution
  + with shared memory (on NUMA platforms)
  + with MPI (for inter- and intra-node execution)
+ gas models
  + ideal gas
  + mixtures of thermally perfect gases
+ kinetics
  + finite-rate chemistry (for thermally perfect gas mixtures)
+ turbulence models
  + Spalart-Allmaras
  + k-$\omega$
+ conjugate heat transfer (coupling fluid/solid domains)
  + structured grids in 2D or 3D
+ shock-fitting
+ user-defined run-time customisation for:
  + initial conditions
  + boundary conditions	
  + source terms
  + grid motion
+ block-marching mode








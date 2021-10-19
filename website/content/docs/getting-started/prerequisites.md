---
title: "Prerequisite software"
description: "Prerequisite software for gdtk"
lead: ""
date: 2021-05-26
lastmod: 2021-05-26
draft: false
images: []
menu:
  docs:
    parent: "getting-started"
weight: 10
toc: false
---


Our main development environment is Linux but the programs can be deployed on
Linux, flavours of Unix such as MacOS-X, and MS-Windows using WSL2.

The core Eilmer, and L1d solvers and their modules are mainly written in the
D programming language for speed and the benefits of compile-time checking.
The pre- and post-processing modes make use of the Lua scripting language
so that we get flexibility and convenient customization.
There is also some Ruby a little Tcl/Tk used in the automated testing scripts.

To run simulations, you need an executable versions of the Eilmer and/or L1d programs.
You may build these executable programs from the source code, as described below.


## Prerequisites: building from source

The main requirement is a D language compiler.
We recommend using the latest stable release of the LLVM D compiler.
*Do not install the LLVM D compiler that is bundled in your package manager.*
These are typically quite stale and will fail to build the code.
Instead, follow our [install-by-hand guide.]({{<relref install-d-compiler>}})

To build Eilmer and other programs in the toolkit, you will require:

  + D compiler
      + Binary releases for the latest stable release of the LLVM D compiler (`ldc2` and `ldmd2`)
        may be found at: <https://github.com/ldc-developers/ldc/releases> .
        An install guide for the LLVM D compiler is available [here]({{< relref install-d-compiler >}}).
  + A C compiler
      + GNU compiler is a good option and comes standard on most systems.
  + The gfortran compiler (and 32-bit libraries)
      + gfortran and gfortran-multilib on Debian/Ubuntu/Mint
      + gcc-gfortran on RedHat/CentOS/Fedora
  + git (to clone the repository)
  + readline development package:
      + libreadline-dev on Debian/Ubuntu/Mint
      + readline-devel on RedHat/CentOS/Fedora
  + ncurses development package:
      + libncurses5-dev on Debian/Ubuntu/Mint
      + ncurses-devel on RedHat/CentOS/Fedora
  + openmpi development package:
      + libopenmpi-dev on Debian/Ubuntu/Mint
      + openmpi-devel on RedHat/CentOS/Fedora
        (after install on RedHat-family systems, load with `module load mpi/openmpi-x86_64`,
        and you might like to place that in your `.bashrc` file so that it's loaded every
        time you start a session)
  + plotutils development package:
      + libplot-dev on Debian/Ubuntu/Mint
      + plotutils-devel on RedHat/CentOS/Fedora (for CentOS 8.x, enable PowerTools repo)
  + foreign-function interface packages for Python and Ruby:
      + python3-cffi on Debian/Ubuntu/Mint and RedHat/CentOS/Fedora
      + ruby-ffi on Debian/Ubuntu/Mint

Additionally, if you want to run the test suite, you will require:

  + Ruby package
  + TCL package
  + the Python sympy package

For viewing and plotting results, we recommend:

  + Paraview
  + Gnuplot

The source code of the Lua interpreter is included in the source code repository.


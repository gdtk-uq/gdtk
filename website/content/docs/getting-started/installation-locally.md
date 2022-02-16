---
title: "Install on linux laptop/desktop/VM"
description: "Install on linux laptop/desktop/VM"
lead: ""
date: 2021-05-26
lastmod: 2021-05-26
draft: false
images: []
menu:
  docs:
    parent: "getting-started"
weight: 20
toc: true
---

> **Note:** Be sure that you have installed the prerequisite software before
> proceeding with these instructions.

## Getting the source code

The full source code for the toolkit programs, including a set of examples,
can be found in a public repository on [github](https://github.com/gdtk-uq/gdtk).
To get your own copy, use the git revision control client to clone the repository
with something like the following command:

    git clone https://github.com/gdtk-uq/gdtk.git gdtk

and within a couple of minutes, depending on the speed of your network connection,
you should have your own copy of the full source tree and the complete repository history.


## Installing Eilmer

The default installation directory is `$HOME/gdtkinst`.
To compile and install Eilmer, move into the eilmer source
area and use `make` to coordinate the compiling and installing:

    cd gdtk/src/eilmer
    make install

If you are on a Mac, you'll need to give the `make` command an
extra hint:

    make PLATFORM=macosx install


## Installing the loadable library (optional)

The loadable library consists of the gas module functions compiled into
a dynamically loadable library and a collection of Python and Ruby modules
for loading that library.
Once installed, you should be able to use the gas functions from within
a Python or Ruby interpreter.

To compile and install the library and its supporting wrappers,
move to the gas source directory and use `make` again.

    cd gdtk/src/gas
    make install

Note that the loadable library needs to be built with the DMD64 compiler
and that you need the Foreign-Function-Interface extensions for your
Python and Ruby interpreters.
On a LinuxMint system these packages are `python-cffi` and `ruby-ffi`.


## Installing L1d (optional)

L1d installs into the same location as Eilmer but, its source code
is in a different location.
To compile and install L1d, move into its source area and, again,
use `make` to coordinate the compiling and installing:

    cd gdtk/src/l1d
    make install


## Setting environment variables

We'll assume you are happy using the default install area `$HOME/gdtkinst`.
The next step is to configure your environment to use Eilmer.
You will need to set the
`DGD` variable to point to the top of the installation tree, and
the `DGD_REPO` variable to point to the top of the repository tree.
Note that the installation tree and repository tree are separate. You then
also need to set `$PATH`, `$DGD_LUA_PATH` and `$DGD_LUA_CPATH` to point to
the appropriate places. Some example lines from a .bashrc file are:

    export DGD=$HOME/gdtkinst
    export DGD_REPO=$HOME/gdtk
    export PATH=$PATH:$DGD/bin
    export DGD_LUA_PATH=$DGD/lib/?.lua
    export DGD_LUA_CPATH=$DGD/lib/?.so


To use the gas models as a loadable library in Python and Ruby,
you should set the following environment variables, also.

    export PYTHONPATH=${PYTHONPATH}:${DGD}/lib
    export RUBYPATH=${RUBYPATH}:${DGD}/lib
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DGD}/lib


Remember to refresh your current shell (or log out and log back in) so
that your newly configured environment is available.






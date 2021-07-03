---
title: "Install on a Mac"
date: 2018-08-06T08:05:42+10:00
draft: true
---

These notes were prepared by Ingo Jahn.
They document the procedure Ingo used to prepare
a development environment for building and running
the codes in the DGD collection.

## Steps to prepare Mac for Eilmer use

1. Install Xcode developer tools. Open terminal and type:

    $ xcode-select -p

2. Install newest version of python (optional)

    https://www.python.org/downloads/mac-osx/

3. Install homebrew

    https://brew.sh/

4. Install gedit editor

    $ brew install gedit

5. Install git

    $ brew install git

6. Install dmd

    $ brew install dmd

7. Create `dgd` directory and clone repository.  See [Getting the source code](/docs/getting-started/#getting-the-source-code).

8. Install a copy of gfortran (included as part of gcc)

    $ brew install gcc

9. (For MPI use) Install OpenMPI v 1.8 from source.

    This is an important note. The version of MPI bundled with OSX is too new to
    work with the DLang openmpi bindings as used by Eilmer. So if you use the
    OSX-bundled OpenMPI, you will not be able to build Eilmer. Do not install
    openmpi from the Max packages.

    Instead, you can download and install an old version of openmpi. We have
    tested 1.8 and that works well. It's not uncommon to have to build your
    own version of openmpi on Mac. There are more instructions available
    at this openmpi FAQ link:

    https://www.open-mpi.org/faq/?category=osx#osx-bundled-ompi

10. Build Eilmer by following instructions [here](/docs/getting-started). NOTE: PLATFORM=osx 

11. Install Paraview
   
    https://www.paraview.org/download/



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
toc: true
---


Our main development environment is Linux but the programs can be deployed on
Linux, flavours of Unix such as MacOS, and MS-Windows using WSL2.

The core Eilmer, and L1d solvers and their modules are mainly written in the
D programming language for speed and the benefits of compile-time checking.
The pre- and post-processing modes make use of the Lua scripting language
so that we get flexibility and convenient customization.
There is also some Ruby, Python and a little Tcl/Tk used in the automated testing scripts.

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
  + Make and a C compiler
      + GNU compiler is a good option and comes standard on most systems.
      + On Debian/Ubuntu/Mint, install the package `build-essential`.
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



## Notes on Windows 10 and WSL 2

Before installing Eilmer, you need to get a Linux environment set up on your Windows 10 computer.

### Linux installation on WSL2

If you have a sufficiently recent version of Windows 10,
follow these instructions to install a Linux distribution:
[Windows Subsystem for Linux Installation Guide for Windows 10](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
This approach uses the Microsoft Store to get the Linux distribution.

Instead of going through the Microsoft Store, you may want to down a WSL Linux distribution manually.
[These notes](https://docs.microsoft.com/en-us/windows/wsl/install-manual) will explain how to do so.

If you are new to Linux, you could follow [WSL](https://wiki.ubuntu.com/WSL) notes
from the Ubuntu world to get a bit of an introduction to running Linux.
The notes describe typical tasks such as installing packages in Ubuntu Linux.

### Mounting Windows drives

From [this blog entry](https://docs.microsoft.com/en-us/archive/blogs/wsl/file-system-improvements-to-the-windows-subsystem-for-linux)
we get the following example of mounting an external drive using DrvFs.

To mount a removable drive D: as /mnt/d directory, run the following commands:

    sudo mkdir /mnt/d
    sudo mount -t drvfs D: /mnt/d

Now, you will be able to access the files of your D: drive under /mnt/d.
When you wish to unmount the drive, for example so you can safely remove it, run the following command:

    sudo umount /mnt/d


### Xserver for Windows (optional)

To use graphical applications install an Xserver such as [VcXsrv](https://sourceforge.net/projects/vcxsrv/).
There are some Windows firewall settings that you will need to adjust and
you will need to have an appropriate value for the DISPLAY variable in your Linux environment.
See the notes on getting the Xserver working from WSL2:
[Running an X server with WSL2](https://skeptric.com/wsl2-xserver/)

Most of our programs, such as Eilmer and L1d, are console-based
so you do not need an Xserver running to use them.
GUI-based programs, such as Paraview, will need an Xserver to run.

### Install Eilmer On WSL2

With a Linux distribution now running on your Windows workstation,
you now can continue with `Getting the source code` in
the [getting-started]({{< relref "quick-start" >}}) instructions.


## Notes on MacOS

Eilmer runs nicely on a Mac Mini with the Apple M1 chip, using MacOS Big Sur.

To get the necessary compilers and libraries onto your Mac, start by opening a terminal
and installing the Command Line Tools for Xcode.

    xcode-select --install

and then install the Homebrew package manager by following the installation directions
at [https://brew.sh](https://brew.sh).
It's just a copy and paste of a single command into your terminal.

With the Homebrew package manager installed, you can install the ldc compiler and dub from the terminal.
Other Homebrew packages to install are `gcc` (to get `gfortran`), `open-mpi`, `plotutils`, `readline`, `ncurses`,
`python`, `tcl-tk` and `gnuplot`.
The package manager will determine if other dependencies need to be installed.

    brew install ldc
    brew install dub
    brew install gcc
    brew install open-mpi
    brew install readline
    brew install ncurses
    brew install python
    brew install tcl-tk
    brew install gnuplot
    brew install --cask paraview

Note that the newly installed gcc will be something like gcc-11,
which is much newer than the system-installed gcc on Big Sur.
Later, during the Eilmer build process, we see the linker complain about compiler modules being
built by older versions of the compiler but it all seems to go together and work.

Also, the Homebrew Python (invoked from the command line as `python3`)
will be newer than the system-installed python so you should install its extension packages
using `pip3`.

    pip3 install numpy
    pip3 install matplotlib
    pip3 install sympy
    pip3 install cffi

The command interpreter in the terminal is the Z shell so, when setting environment variables,
put your `export` statements into `.zshrc`
You can now continue with `Getting the source code` in the
[getting-started]({{< relref "quick-start" >}}) instructions.



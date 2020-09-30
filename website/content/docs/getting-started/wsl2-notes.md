---
title: "WSL2 notes"
date: 2020-09-30
weight: 40
draft: false
---

## Notes on Windows 10 and WSL 2

This page provides some pointers to documents on
how to get a Linux environment set up on your Windows 10 computer.

### Linux installation

Follow these instructions to install a Linux distribution:
[Windows Subsystem for Linux Installation Guide for Windows 10](https://docs.microsoft.com/en-us/windows/wsl/install-win10)

If you are new to Linux, you could follow [WSL](https://wiki.ubuntu.com/WSL) notes
from the Ubuntu world to get a bit of an introduction to running Linux.
The notes describe typical tasks such as installing packages in Ubuntu Linux.


### Xserver (optional)

To use graphical applications install an Xserver such as [VcXsrv](https://sourceforge.net/projects/vcxsrv/).
There are some Windows firewall settings that you will need to adjust and
you will need to have an appropriate value for the DISPLAY variable in your Linux environment.
See the notes on getting the Xserver working from WSL2:
[Running an X server with WSL2](https://skeptric.com/wsl2-xserver/)

Most of our programs, such as Eilmer and L1d, are console-based
so you do not need an Xserver running to use them.
GUI-based programs, such as Paraview, will need an Xserver to run.


### DMD compiler

With your favourite Linux distribution now running,
it's time to install the DMD compiler.
You may choose a prebuilt binary for your selected Linux distribution or
you can install the DMD compiler via the install script at https://dlang.org/download.html

    curl -fsS https://dlang.org/install.sh | bash -s dmd

Having installed the compier via the script,
you will need to activate the DMD compiler in order to use it.
At the time of writing these notes, the current DMD compiler version was 2.094.0.
To activate the compiler, issue the command:

    source ~/dlang/dmd-2.094.0/activate

and, to deactivate:

    deactivate


### Install Eilmer

With Linux running on your workstation,
you can continue with the [getting-started]({{< relref "quick-start" >}}) instructions.



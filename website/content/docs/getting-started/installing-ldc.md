---
title: "Guide to Installing LDC"
date: 2020-04-01T08:22:59+10:00
weight: 20
draft: false
---

## Installing the LDC compiler
You need a D compiler to build Eilmer from source and we recommend
the LDC compiler.
If you have some linux chops, the install process is quite straightforward:

  + Download the latest stable release from:
  https://github.com/ldc-developers/ldc/releases
  + Unpack in an area of your choosing.
  + Set your `PATH` variable so that it can find the binaries in the `bin/` subdirectory
    of what you just unpacked.

For those newer to linux, read on for a step-by-step walk through.
Advanced users might just skip to the [bottom step]({{< relref "#for-advanced-users-of-the-loadable-library" >}}) about configuring
your environment to use the loadable library.

### Step-by-step walk through for installing the LDC compiler

For this example, I will install version `ldc-1.20.1`.
You may find a newer stable release when you visit the ldc webpage.
That is fine.
In fact, we often encourage use of the newest stable release.
For most of this you will work in a terminal.
The first step involves the browser to find a package to download.

1. Go to the LDC project github page ( https://github.com/ldc-developers/ldc/releases ),
   scroll down to find the latest stable release.
   (Note: sometimes a pre-release is listed first. Skip this. Go for a stable release.)

2. Download the binary package for x86_64 linux. It will have the form:
   `ldc2-X.YY.Z-linux-x86_64.tar.xz`.
   For this example, I download: ldc2-1.20.1-linux-x86_64.tar.xz.

   This file is what we call a compressed tarball. The tarball refers to the fact
   that all of the files in this package have been grouped together into one file.
   The compressed part means that tarball has been modified with a compression program.
   In this case, the `xz` compression program was used.

   > For the remainder of these steps, open a terminal.

3. Copy the compressed tarball to an area for installation. I like to use
   an `opt/` directory for extra programs that I install by hand in this manner.
   The following are terminal commands you can type:

        cd
        mkdir opt
        cd opt
        cp ~/Downloads/ldc2-1.20.1-linux-x86_64.tar.xz .

4. Now unpack the compressed tarball.

        tar -Jxvf ldc2-1.20.1-linux-x86_64.tar.xz

   You will see it print a list of all the files that have been unpacked.

   This command has actually done two steps in one: decompressed the file
   and extracted the files from the large archive file. The `tar` command
   is responsible for the extraction, and we gave it the hint to do the
   decompression first with the `-J` flag.

   There are in fact four flags we used here: `J`, `x`, `v`, and `f`.
   They were grouped together and passed as one to the command by
   prepending the dash, ie. `-Jxvf`.
   The `J` is used to tell `tar` to outsource the decompression to the `xz` program.
   The `x` tells `tar` to do extraction.
   The `v` asks the command to be a little verbose and tell us what it's doing.
   Without that flag, the command still works, but you would see no output on
   your screen.
   The `f` tells `tar` to work on a file. This might seem odd -- don't we always work
   on files? Well, `tar` stands for tape archive. It was written in the days of
   computers when files were archived onto tape. To this day, `tar`'s default
   action will be to try to put things onto tape, and read from tape.
   So you will very very often see the `-f` flag in usage on modern computers since
   we are often working with files on local disk.

5. Set your `PATH` variable so that it finds the LDC binaries.

   To do this, I add a line to my `.bashrc` file, which is a hidden file
   sitting in my home directory.
   On a Debian-based system, the preferred place for modifications is
   in the file `.bash_aliases`.
   If you are using Ubuntu or Mint Linux, then you are using a Debian-based system.
   You will need to open an editor to make the addition to that file.
   Let's try `gedit` in this example, but any editor will do.

        cd
        gedit .bashrc

   Now add a line at the end of the file. Perhaps you did not download
   version 1.20.1. I did, so I'll use that in this example. Adjust
   your line according to the version you downloaded.

        export PATH=${PATH}:${HOME}/opt/ldc2-1.20.1-linux-x86_64/bin

   The `PATH` variable is what we call an environment variable.
   It configures a particular aspect of your working environment.
   This particular variable tells the compute which directories
   of all the many directories on your hard drive to search
   when looking for binaries.
   This is an efficiency measure.
   If the computer searched every directory every time you typed
   a command, the response time on commands would be very slow.
   Instead, we configure the `PATH` variable so that it only contains
   a very small subset of all the total directories to search.
   This speeds up enormously the search for commands when we type
   them at the command prompt.

   Note that we *appended* the new search directory to what
   was already set in `PATH`.
   We did that with the expression
   `${PATH}:${HOME}/opt/ldc2-1.20.1-linux-x86_64/bin`.
   This is important.
   Using the append expression, we preserve the system default `PATH` search
   directories.
   If we were to override that, we would lose access to basic
   commands and programs on our system, all the things we think
   of as system defaults.

6. Changes to your `.bashrc` file do not take immediate effect in this
   session. To have your `.bashrc` file re-processed, do:

        cd
        . .bashrc

   Look at that second command carefully. It is strange at first glance.
   It is a dot followed by a space followed by the filename `.bashrc`.
   It is an instruction to process what appears in `.bashrc`.

### For advanced users of the loadable library

If you are using the loadable library, then you will also need to configure your
`LD_LIBRARY_PATH`. Again make the appropriate change in `.bashrc` or equivalent:

        export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/opt/ldc2-1.20.1-linux-x86_64/lib


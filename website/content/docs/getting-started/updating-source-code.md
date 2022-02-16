---
title: "Updating the source code"
description: "Updating the source code"
lead: ""
date: 2021-05-26
lastmod: 2021-05-26
draft: false
images: []
menu:
  docs:
    parent: "getting-started"
weight: 50
toc: false
---

Because Eilmer, L1d and the toolkit functions are being actively developed,
we will make frequent changes
and additions to the source code and the examples.
To update your copy of the repository,
move into any directory within the working tree and pull any new revisions into your local copy.

    cd gdtk/src/eilmer
    make clean
    git pull -v

Then, to build a refreshed copy of Eilmer,
use `make` to coordinate the compiling and installing as before:

    make install





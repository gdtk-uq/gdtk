TinyEndian
==========

![](https://travis-ci.org/dlang-community/tinyendian.svg?branch=master) ![](https://img.shields.io/dub/v/tinyendian.svg)

Introduction
------------

TinyEndian is a minimal endianness library for the D programming
language. It has no external dependencies, it only needs a D compiler
and Phobos (standard library). TinyEndian doesn't allocate memory and is
fully `@nogc` to allow use in high-performance code.

The API is not stable and may change in the future.

Features
--------

-   Swap byte order of 2- or 4-byte elements in an array in place.
-   Read a UTF-8, UTF-16 or UTF-32 buffer, determine its endianness
    using a UTF byte-order-mark and convert it to system endianness in
    place.
-   No external dependencies.
-   pure, nothrow and @nogc.

Directory structure
-------------------

| Directory  | Contents                                     |
|------------|----------------------------------------------|
| `./`       | This README file, license, DUB package file. |
| `./source` | Source code.                                 |

Usage
-----

Assuming you use [dub](http://code.dlang.org/about), add this line:

    "tinyendian": { "version" : "~>0.2.0" }

to the `"dependencies"` in your project's `dub.json`.

If you don't use dub, you can directly copy the `source/tinyendian.d`
file into your project.

TinyEndian requires DMD 2.067 or better.

License
-------

TinyEndian is released under the terms of the Boost Software License
1.0. This license allows you to use the source code in your own
projects, open source or proprietary, and to modify it to suit your
needs. However, in source distributions, you have to preserve the
license headers in the source code and the accompanying license file.

Full text of the license can be found in file `LICENSE_1_0.txt` and is
also displayed here:

    Boost Software License - Version 1.0 - August 17th, 2003

    Permission is hereby granted, free of charge, to any person or organization
    obtaining a copy of the software and accompanying documentation covered by
    this license (the "Software") to use, reproduce, display, distribute,
    execute, and transmit the Software, and to prepare derivative works of the
    Software, and to permit third-parties to whom the Software is furnished to
    do so, all subject to the following:

    The copyright notices in the Software and this entire statement, including
    the above license grant, this restriction and the following disclaimer,



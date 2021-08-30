---
title : "foamMesh"
description: "foamMesh"
lead: "grid generation for OpenFOAM"
date: 2021-08-30
lastmod: 2021-08-30
draft: false
images: []
menu:
  docs:
    parent: "grid-utilities"
weight: 10
toc: false
---


**foamMesh** is a stand-alone program and an extension of the geometry package developed
for the Eilmer compressible-flow simulation program.
It allows gas flow and solid domains, generated using the Eilmer geometry
package, to be converted into grids suitable for OpenFOAM simulations.

The generation of OpenFOAM grids is achieved by adding an extra step
to the Lua scripts, used to convert the 2D and 3D grids of finite volume
cells into corresponding foam meshes.
At the same stage labels are assigned to the outward facing edges in the
(x,y)-plane for 2D grids, or outward facing patches for 3D grids, which allow
definition of the OpenFoam boundary conditions.

## User guide
<div class="row">
   <div class="col-sm-5">
   <a href="/pdfs/foammesh-user-guide.pdf">
    <img src="/images/foammesh-thumbnail.png" style="width:100%">
   </a>
   </div>
   <div class="col-sm-10">
   The <a href="/pdfs/foammesh-user-guide.pdf"> foamMesh User Guide</a> is available as PDF.
   </div>
</div>



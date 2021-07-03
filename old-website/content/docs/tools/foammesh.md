---
title: "foamMesh"
draft: false
weight: 50
---

# foamMesh: grid generation for OpenFOAM

**foamMesh** is a stand-alone program and an extension of the geometry package developed
for the Eilmer4 compressible-flow simulation program.
It allows gas flow and solid domains, generated using the Eilmer4 geometry
package, to be converted into grids suitable for OpenFOAM simulations.

The generation of OpenFOAM grids is achieved by adding an extra step
to the Lua scripts, used to convert the 2D and 3D grids of finite volume
cells into corresponding foam meshes.
At the same stage labels are assigned to the outward facing edges in the
(x,y)-plane for 2D grids, or outward facing patches for 3D grids, which allow
definition of the OpenFoam boundary conditions.

## User guide
[User guide in PDF](/pdfs/foammesh-user-guide.pdf)




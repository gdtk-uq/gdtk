// paver2d.d
// An unstructured grid generator for two-dimensional regions.
// Based on Heather Muir's Paver that she built for her UG thesis in 2015,2016.
// PJ, 2018-02-26
//

module paver2d;

import std.stdio;
import geom;

void fill_interior(ref Vector3[] vertices, ref USGFace[] faces, ref USGCell[] cells)
{
    // Generate quadrilateral cells for a closed region.
    // The discretized boundary is already stored in the vertices array and
    // the faces on the boundary are already present in the faces array.
    // Ther are no cells at this point.
    writeln("[TODO] Finish code");
} // end fill_interior()

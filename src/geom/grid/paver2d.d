// paver2d.d
// An unstructured grid generator for two-dimensional regions.
// Based on Heather Muir's Paver that she built for her UG thesis in 2015,2016.
// PJ, 2018-02-26
//

module geom.grid.paver2d;

import std.stdio;
import geom;

enum NodeType {side, rowend, corner, reversal}

void fill_interior(ref Vector3[] vertices,
                   ref USGFace[] faces, ref size_t[string] faceIndices,
                   ref USGCell[] cells,
                   ref size_t[] bndry_vtx_ids)
{
    writeln("Cheat by making all triangle cells connected to the mid-point.");
    // This is just a temporary arrangement, so that we can generate
    // any sort of grid while we get to understand Heather's code.
    size_t n = bndry_vtx_ids.length;
    Vector3 mid = Vector3(0,0);
    foreach (id; bndry_vtx_ids) { mid += vertices[id]; }
    mid /= n;
    size_t id_mid = vertices.length;
    vertices ~= mid;
    // Make all of the interior faces.
    foreach (id; bndry_vtx_ids) {
        size_t[] vtx_id_list = [id, id_mid];
        string faceTag = makeFaceTag(vtx_id_list);
        faceIndices[faceTag] = faces.length;
        faces ~= new USGFace(vtx_id_list);
    }
    // Make the actual cells as triangles with all faces pointing out.
    foreach (j; 0 .. n-1) {
        size_t id_j = bndry_vtx_ids[j];
        size_t id_jp1 = (j+1 < n) ? bndry_vtx_ids[j+1] : bndry_vtx_ids[0];
        size_t[] vtx_id_list = [id_j, id_jp1, id_mid];
        size_t[] face_id_list = [faceIndices[makeFaceTag([id_j, id_jp1])],
                                 faceIndices[makeFaceTag([id_jp1, id_mid])],
                                 faceIndices[makeFaceTag([id_mid, id_j])]];
        size_t id = cells.length;
        cells ~= new USGCell(USGCell_type.triangle, id, vtx_id_list, face_id_list, [1,1,1]);
    }
    
    writeln("[TODO] Write the actual paver code");
    // Generate quadrilateral cells for a closed region.
    // The discretized boundary is already stored in the vertices array and
    // the faces on the boundary are already present in the faces array.
    // Ther are no cells at this point.
} // end fill_interior()

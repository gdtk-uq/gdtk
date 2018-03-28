/** grid_deform.d
 * This file holds the methods used to deform a mesh for shape optimisation.
 *
 * Author: Kyle D.
 * Date: 2018-01-17
**/

module grid_deform;

import std.stdio;
import fluidblock;
import std.math;
import std.algorithm;
import geom;

void inverse_distance_weighting(FluidBlock blk, size_t gtl, string[] NonFixedBoundaryList=[]) {
    
    // we assume the boundary vertices have already been perturbed, and have their new positions stored at the current gtl.
    // now we compute perturbed internal vertices
    foreach(vtxi; blk.vertices) {
	double[2] numer_sum; numer_sum[0] = 0.0; numer_sum[1] = 0.0; double denom_sum = 0.0; 
	foreach(bndary; blk.bc) {
            if ( NonFixedBoundaryList.canFind(bndary.group) == false) {
                foreach(vtxs; bndary.vertices) {
                    Vector3 ds;
                    ds.refx = vtxs.pos[gtl].x - vtxs.pos[0].x;
                    ds.refy = vtxs.pos[gtl].y - vtxs.pos[0].y;
                    double r; double w;
                    double dx = vtxi.pos[0].x - vtxs.pos[0].x; 
                    double dy = vtxi.pos[0].y - vtxs.pos[0].y; 
                    r = sqrt(dx*dx + dy*dy);
                    w = 1.0 / (r*r);
                    numer_sum[0] += ds.x * w; numer_sum[1] += ds.y * w;
                    denom_sum += w;
                }
            }
        }
	// update vertex positions
	if (blk.boundaryVtxIndexList.canFind(vtxi.id) == false) {
	    vtxi.pos[gtl].refx = vtxi.pos[0].x + (numer_sum[0]/denom_sum);
	    vtxi.pos[gtl].refy = vtxi.pos[0].y + (numer_sum[1]/denom_sum);
	}
    }
}

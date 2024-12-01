// face.d -- Part of the Puffin and Lorikeet flow calculators.
//
// PA Jacobs
// 2022-01-22
//
module face;

import std.format;
import std.math;

import geom;
import gas;
import gasdyn.gasflow;
import config;
import flow;
import cell;


class Face2D {
public:
    CQIndex cqi;
    Vector3 pos;
    Vector3* p0, p1; // pointers to vertices at each end of face
    Vector3 n; // unit normal (to right when looking from p0 to p1)
    Vector3 t1; // unit tangent (from p0 to p1)
    double area; // per unit depth for 2D planar, per radian for axisymmetric
    //
    double[] F; // Flux vector
    //
    Cell2D[2] left_cells; // References to cells on the left, starting with closest.
    Cell2D[2] right_cells; // References to cells on the right, starting with closest.
    //
    // Workspace for the Osher-type flux calculator.
    GasState stateLstar, stateRstar, stateX0;

    this(GasModel gmodel, in CQIndex cqi)
    {
        this.cqi = CQIndex(cqi);
        pos = Vector3();
        F.length = cqi.n;
        // Workspace for Osher-type flux calculator.
        stateLstar = GasState(gmodel);
        stateRstar = GasState(gmodel);
        stateX0 = GasState(gmodel);
    }

    this(ref const(Face2D) other)
    {
        cqi = CQIndex(other.cqi);
        pos = Vector3(other.pos);
        F = other.F.dup;
        stateLstar = GasState(other.stateLstar);
        stateRstar = GasState(other.stateRstar);
        stateX0 = GasState(other.stateX0);
    }

    override
    string toString() const
    {
        string repr = "Face2D(";
        repr ~= format("p0=%s", ((p0) ? (*p0).toString() : "null"));
        repr ~= format(", p1=%s", ((p1) ? (*p1).toString() : "null"));
        repr ~= format(", pos=%s, n=%s, t1=%s, area=%g", pos, n, t1, area);
        repr ~= format(", F=%s", F);
        repr ~= ")";
        return repr;
    }

    @nogc
    void compute_geometry(bool axiFlag)
    // Update the geometric properties from vertex data.
    {
        t1 = *p1; t1 -= *p0; t1.normalize();
        Vector3 t2 = Vector3(0.0, 0.0, 1.0);
        cross(n, t1, t2); n.normalize();
        // Assume unit depth in z for 2D planar geometry.
        area = distance_between(*p1, *p0);
        pos = *p0; pos += *p1; pos *= 0.5;
        // Area per radian about the z-axis for 2D axisymmetric geometry.
        if (axiFlag) { area *= pos.y; }
        return;
    }

} // end class Face2D

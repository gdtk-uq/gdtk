// cell.d -- Part of the Puffin steady-flow calculator.
//
// PA Jacobs
// 2022-01-22
//
module cell;

import std.format;

import geom;
import gas;
import config;
import flow;

enum Compass {north=0, east=1,    south=3, west=4};
// alias      upper    downstream lower    upstream

class Cell2D {
public:
    Vector3 pos;
    Vector3* p00, p10, p11, p01;
    Face[] faces;
    double volume, xyplane_area;
    double iLen, jLen, kLen; // distances across the cell
    FlowState2D fs;

    this(GasModel gm)
    {
        pos = Vector3();
        fs = new FlowState2D(gm);
    }

    this(ref const(Cell2D) other)
    {
        pos = Vector3(other.pos);
        fs = new FlowState2D(other.fs);
    }

    override
    string toString() const
    {
        string repr = "Cell2D(";
        repr ~= format("pos=%s", pos);
        repr ~= format(", fs=%s", fs);
        repr ~= ")";
        return repr;
    }

    @nogc
    int compute_geometry()
    // Update the geometric properties from vertex data.
    {
        xyplane_quad_cell_properties(*p00, *p10, *p11, *p01, pos, xyplane_area, iLen, jLen, kLen);
        volume = xyplane_area * ((Config.axisymmetric) ? pos.y : 1.0);
        return 0;
    }
} // end class Cell2D


class Face2D {
public:
    Vector3 pos;
    Vector3* p0, p1; // pointers to vertices at each end of face
    Vector3 n; // unit normal (to right when looking from p0 to p1)
    Vector3 t1; // unit tangent (from p0 to p1)
    double area; // per unit depth for 2D planar, per radian for axisymmetric

    this()
    {
        pos = Vector3();
    }

    this(ref const(Face2D) other)
    {
        pos = Vector3(other.pos);
    }

    override
    string toString() const
    {
        string repr = "Face2D(";
        repr ~= format("pos=%s", pos);
        repr ~= ")";
        return repr;
    }

    @nogc
    int compute_geometry()
    // Update the geometric properties from vertex data.
    {
        n = *p1; n -= *p0; n.normalize();
        Vector3 t2 = Vector3(0.0, 0.0, 1.0);
        cross(t1, n, t2);
        area = distance_between(*p1, *p0);
        pos = *p0; pos += *p1; pos *= 0.5;
        if (Config.axisymmetric) { area *= pos.y; }
        return 0;
    }
} // end class Face2D

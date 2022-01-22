// cell.d -- Part of the Puffin steady-flow calculator.
//
// PA Jacobs
// 2022-01-22
//
module cell;

import std.format;

import geom;
import gas;
import flow;


class Cell {
public:
    Vector3 pos;
    FlowState2D fs;

    this(GasModel gm)
    {
        pos = Vector3();
        fs = new FlowState2D(gm);
    }

    this(ref const(Cell) other)
    {
        pos = Vector3(other.pos);
        fs = new FlowState2D(other.fs);
    }

    override
    string toString() const
    {
        string repr = "Cell(";
        repr ~= format("pos=%s", pos);
        repr ~= format(", fs=%s", fs);
        repr ~= ")";
        return repr;
    }
} // end class Cell


class Face {
public:
    Vector3 pos;

    this()
    {
        pos = Vector3();
    }

    this(ref const(Face) other)
    {
        pos = Vector3(other.pos);
    }

    override
    string toString() const
    {
        string repr = "Face(";
        repr ~= format("pos=%s", pos);
        repr ~= ")";
        return repr;
    }
} // end class Face


class Vertex {
public:
    Vector3 pos;

    this()
    {
        pos = Vector3();
    }

    this(ref const(Vertex) other)
    {
        pos = Vector3(other.pos);
    }

    override
    string toString() const
    {
        string repr = "Face(";
        repr ~= format("pos=%s", pos);
        repr ~= ")";
        return repr;
    }
} // end class Vertex

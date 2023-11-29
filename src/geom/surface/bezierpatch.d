/*
 * Authors: Rowan G. and Peter J.
 * History: Ported from Eilmer3 on 2018-10-31
 */

module geom.surface.bezierpatch;

import std.stdio : File;
import std.conv;
import ntypes.complex;
import nm.number;

import geom.elements;
import geom.gpath;
import geom.surface.parametricsurface;


class BezierPatch : ParametricSurface {
public:
    // Control net is defined by points Q[i][j] with
    // 0 <= i <= n, 0 <= j <= m
    Vector3[][] Q;
    int n; // order of polynomial in r-parameter direction
    int m; // order of polynomial in s-parameter direction

    this(const Vector3[][] _Q, int _n, int _m)
    {
        n = _n;
        m = _m;
        Q.length = _n + 1;
        foreach (i; 0 .. Q.length) {
            Q[i] = _Q[i].dup;
        }
        initialiseWorkingSpace();
    }

    this(ref const(BezierPatch) other)
    {
        n = other.n;
        m = other.n;
        Q.length = other.Q.length;
        foreach (i; 0 .. Q.length) {
            Q[i] = other.Q[i].dup;
        }
        initialiseWorkingSpace();
    }

    override BezierPatch dup() const
    {
        return new BezierPatch(this.Q, this.n, this.m);
    }

    override Vector3 opCall(double r, double s) const
    {
        Vector3[] B;
        B.length = n + 1;
        foreach (i; 0 .. n + 1) {
            auto p = _b1s[i](s);
            B[i].set(p);
        }
        Bezier b2 = new Bezier(B);
        return b2(r);
    }

    override string toString() const
    {
        return "BezierPatch()";
    }

    void initialiseWorkingSpace()
    {
        _b1s.length = n + 1;
        foreach (i; 0 .. n + 1) {
            _b1s[i] = new Bezier(Q[i]);
        }
    }

private:
    Bezier[] _b1s;
} // end class BezierPatch

void writeCtrlPtsAsVtkXml(BezierPatch bPatch, string fileName)
{
    auto n = bPatch.Q.length;
    auto m = bPatch.Q[0].length;

    auto f = File(fileName, "w");
    f.writeln("<VTKFile type=\"StructuredGrid\" version=\"1.0\" header_type=\"UInt64\">");
    f.writefln("  <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">", 0, n-1, 0, m-1);
    f.writefln("    <Piece Extent=\"%d %d %d %d 0 0\">",  0, n-1, 0, m-1);
    f.writeln("      <Points>");
    f.writeln("        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">");
    foreach (j; 0 .. m) {
        foreach (i; 0 .. n) {
            auto p = bPatch.Q[i][j];
            f.writefln("       %20.16e %20.16e %20.16e", p.x, p.y, p.z);
        }
    }
    f.writeln("        </DataArray>");
    f.writeln("      </Points>");
    f.writeln("    </Piece>");
    f.writeln("  </StructuredGrid>");
    f.writeln("</VTKFile>");
    f.close();
}

version(bezierpatch_test) {
    import util.msg_service;
    int main() {
        double L = 2.0;
        double H = 1.0;
        int n = 3;
        int m = 4;
        double dx = L/n;
        double dy = H/m;
        Vector3[][] Q;
        Q.length = n + 1;
        foreach (i; 0 .. n + 1) {
            Q[i].length = m + 1;
            foreach (j; 0 .. m + 1) {
                Q[i][j] = Vector3(i*dx, j*dy, 0.0);
            }
        }
        auto bezPatch = new BezierPatch(Q, n, m);
        auto p = bezPatch(0.2, 0.25);
        assert(approxEqualVectors(p, Vector3(0.4, 0.25, 0.0)), failedUnitTest());
        p = bezPatch(0.9, 0.1);
        assert(approxEqualVectors(p, Vector3(1.8, 0.1, 0.0)), failedUnitTest());
        return 0;
    }
}

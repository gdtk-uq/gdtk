// tfivolume.d

module geom.volume.tfivolume;

import std.conv;
import std.format;
import ntypes.complex;
import nm.number;

import geom.geometry_exception;
import geom.elements;
import geom.gpath;
import geom.surface;
import geom.volume.parametricvolume;


class TFIVolume : ParametricVolume {
public:
    ParametricSurface[6] faces;
    Vector3[8] p;
    // Path[12] edges;

    // Generic, 6-face constructor
    this(in ParametricSurface[] faceArray)
    {
        // Keep a local copy of the defining faces and parametric ranges.
        foreach(i; 0 .. 6) faces[i] = faceArray[i].dup();
        // Check that the corners of the faces coincide.
        p[0] = faces[Face.bottom](0.0,0.0);
        Vector3 p_alt1 = faces[Face.west](0.0,0.0);
        Vector3 p_alt2 = faces[Face.south](0.0,0.0);
        if (!approxEqualVectors(p[0], p_alt1) || !approxEqualVectors(p[0], p_alt2)) {
            throw new Error(text("TFIVolume open at corner 0 p= ", p[0],
                                 " p_alt1= ", p_alt1," p_alt2= ", p_alt2));
        }
        p[1] = faces[Face.bottom](1.0,0.0);
        p_alt1 = faces[Face.east](0.0,0.0);
        p_alt2 = faces[Face.south](1.0,0.0);
        if (!approxEqualVectors(p[1], p_alt1) || !approxEqualVectors(p[1], p_alt2)) {
            throw new Error(text("TFIVolume open at corner 1 p= ", p[1],
                                 " p_alt1= ", p_alt1," p_alt2= ", p_alt2));
        }
        p[2] = faces[Face.bottom](1.0,1.0);
        p_alt1 = faces[Face.east](1.0,0.0);
        p_alt2 = faces[Face.north](1.0,0.0);
        if (!approxEqualVectors(p[2], p_alt1) || !approxEqualVectors(p[2], p_alt2)) {
            throw new Error(text("TFIVolume open at corner 2 p= ", p,
                                 " p_alt1= ", p_alt1," p_alt2= ", p_alt2));
        }
        p[3] = faces[Face.bottom](0.0,1.0);
        p_alt1 = faces[Face.west](1.0,0.0);
        p_alt2 = faces[Face.north](0.0,0.0);
        if (!approxEqualVectors(p[3], p_alt1) || !approxEqualVectors(p[3], p_alt2)) {
            throw new Error(text("TFIVolume open at corner 3 p= ", p,
                                 " p_alt1= ", p_alt1," p_alt2= ", p_alt2));
        }
        p[4] = faces[Face.top](0.0,0.0);
        p_alt1 = faces[Face.west](0.0,1.0);
        p_alt2 = faces[Face.south](0.0,1.0);
        if (!approxEqualVectors(p[4], p_alt1) || !approxEqualVectors(p[4], p_alt2)) {
            throw new Error(text("TFIVolume open at corner 4 p= ", p[4],
                                 " p_alt1= ", p_alt1," p_alt2= ", p_alt2));
        }
        p[5] = faces[Face.top](1.0,0.0);
        p_alt1 = faces[Face.east](0.0,1.0);
        p_alt2 = faces[Face.south](1.0,1.0);
        if (!approxEqualVectors(p[5], p_alt1) || !approxEqualVectors(p[5], p_alt2)) {
            throw new Error(text("TFIVolume open at corner 5 p= ", p[5],
                                 " p_alt1= ", p_alt1," p_alt2= ", p_alt2));
        }
        p[6] = faces[Face.top](1.0,1.0);
        p_alt1 = faces[Face.east](1.0,1.0);
        p_alt2 = faces[Face.north](1.0,1.0);
        if (!approxEqualVectors(p[6], p_alt1) || !approxEqualVectors(p[6], p_alt2)) {
            throw new Error(text("TFIVolume open at corner 6 p= ", p[6],
                                 " p_alt1= ", p_alt1," p_alt2= ", p_alt2));
        }
        p[7] = faces[Face.top](0.0,1.0);
        p_alt1 = faces[Face.west](1.0,1.0);
        p_alt2 = faces[Face.north](0.0,1.0);
        if (!approxEqualVectors(p[7], p_alt1) || !approxEqualVectors(p[7], p_alt2)) {
            throw new Error(text("TFIVolume open at corner 7 p= ", p[7],
                                 " p_alt1= ", p_alt1," p_alt2= ", p_alt2));
        }
        checkVolume();
    } // end generic constructor

    // Wire-Frame constructor TODO

    // Simple-Box constructor
    this(in Vector3[] p)
    {
        foreach(i; 0 .. 8) { this.p[i] = p[i].dup(); }
        faces[Face.north] = new CoonsPatch(p[3], p[2], p[6], p[7]);
        faces[Face.south] = new CoonsPatch(p[0], p[1], p[5], p[4]);
        faces[Face.east] = new CoonsPatch(p[1], p[2], p[6], p[5]);
        faces[Face.west] = new CoonsPatch(p[0], p[3], p[7], p[4]);
        faces[Face.top] = new CoonsPatch(p[4], p[5], p[6], p[7]);
        faces[Face.bottom] = new CoonsPatch(p[0], p[1], p[2], p[3]);
        // Don't bother testing for open corners.
        checkVolume();
    } // end constructor

    // Surface-extrusion constructor TODO

    // Copy constructor
    this(ref const(TFIVolume) other)
    {
        foreach(i; 0 .. 6) { this.faces[i] = other.faces[i].dup(); }
        foreach(i; 0 .. 8) { this.p[i] = other.p[i].dup(); }
    } // end copy constructor

    override TFIVolume dup() const
    {
        return new TFIVolume(this.faces);
    }

    override Vector3 opCall(double r, double s, double t) const
    // Locate a point within the volume by blended linear interpolation.
    // Input:
    //     r: interpolation parameter i-direction west-->east, 0.0<=r<=1.0
    //     s: interpolation parameter j-direction south-->north, 0.0<=s<=1.0
    //     t: interpolation parameter k-direction bottom-->top, 0.0<=t<=1.0
    // Returns:
    //     a Vector3 value for the point.
    {
        Vector3 pW = faces[Face.west](s,t);
        Vector3 pE = faces[Face.east](s,t);
        Vector3 pS = faces[Face.south](r,t);
        Vector3 pN = faces[Face.north](r,t);
        Vector3 pB = faces[Face.bottom](r,s);
        Vector3 pT = faces[Face.top](r,s);
        double omr = 1.0 - r; double oms = 1.0 - s; double omt = 1.0 - t;
        Vector3 BigC = (omr*oms*omt)*p[0] + (omr*oms*t)*p[4] + (omr*s*omt)*p[3] + (omr*s*t)*p[7] +
            (r*oms*omt)*p[1] + (r*oms*t)*p[5] + (r*s*omt)*p[2] + (r*s*t)*p[6];
        Vector3 p_rst = 0.5*(omr*pW + r*pE + oms*pS + s*pN + omt*pB + t*pT) - 0.5*BigC;
        return p_rst;
    } // end opCall

    override string toString() const
    {
        string repr = "TFIVolume(faces=[" ~ to!string(faces[0]);
        foreach(i; 1 .. 6) { repr ~= ", " ~ to!string(faces[i]); }
        repr ~= "])";
        return repr;
    } // end toString

private:
    void checkVolume()
    {
        // Compute the volume of the hexahedron, as defined by the corners of the block.
        // This should give a rough and ready test for an ill-defined block.
        number volume, iLength, jLength, kLength;
        Vector3 centroid;
        hex_cell_properties(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                            true, centroid, volume, iLength, jLength, kLength);
        if (volume <= 0.0) {
            string msg = format("Apparently-invalid volume (%g) for hex block.", volume);
            throw new GeometryException(msg);
        }
    }
} // end class TFIVolume

version(tfivolume_test) {
    import util.msg_service;
    int main() {
        Vector3[8] p;
        p[0] = Vector3(0.0, 0.1, 0.0);
        p[1] = Vector3(1.0, 0.1, 0.0);
        p[2] = Vector3(1.0, 1.1, 0.0);
        p[3] = Vector3(0.0, 1.1, 0.0);

        p[4] = Vector3(0.0, 0.1, 3.0);
        p[5] = Vector3(1.0, 0.1, 3.0);
        p[6] = Vector3(1.0, 1.1, 3.0);
        p[7] = Vector3(0.0, 1.1, 3.0);

        auto simple_box = new TFIVolume(p);
        auto c = simple_box(0.1, 0.1, 0.5);
        assert(approxEqualVectors(c, Vector3(0.1, 0.2, 1.5)), failedUnitTest());

        ParametricSurface[6] my_faces;
        my_faces[Face.north] = new CoonsPatch(p[3], p[2], p[6], p[7]);
        my_faces[Face.south] = new CoonsPatch(p[0], p[1], p[5], p[4]);
        my_faces[Face.east] = new CoonsPatch(p[1], p[2], p[6], p[5]);
        my_faces[Face.west] = new CoonsPatch(p[0], p[3], p[7], p[4]);
        my_faces[Face.top] = new CoonsPatch(p[4], p[5], p[6], p[7]);
        my_faces[Face.bottom] = new CoonsPatch(p[0], p[1], p[2], p[3]);
        auto generic_box = new TFIVolume(my_faces);
        auto d = generic_box(0.1, 0.1, 0.5);
        assert(approxEqualVectors(d, Vector3(0.1, 0.2, 1.5)), failedUnitTest());
        return 0;
    }
} // end tfivolume_test

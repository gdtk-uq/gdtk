// spherepatch.d
// Generate ParametricSurface patches on a sphere.
// based on the formula for mapping from a cube to a sphere:
//
// Authors: Peter J. and Ingo Jahn
// 2021-04-01 Built in D by adapting Ingo's Lua script.

module geom.surface.spherepatch;

import std.conv;
import std.algorithm;
import std.string;
import std.math;
import nm.complex;
import nm.number;

import geom.elements;
import geom.gpath;
import geom.surface.parametricsurface;


class SpherePatch : ParametricSurface {
public:
    number radius;
    Vector3 centre;
    string face_name; // face identifier as per the canonical cube
    // If we want one half or the other of a full-face patch,
    // we specify the adjoining face nearest the half we want.
    string which_half;

    this(number radius, Vector3 centre, string face_name, string which_half="")
    {
        this.radius = radius;
        this.centre = Vector3(centre);
        auto valid_names = ["north", "east", "south", "west", "bottom", "top"];
        if (!canFind(valid_names, face_name)) {
            throw new Exception("Invalid face name: " ~ face_name);
        }
        this.face_name = face_name;
        if (which_half.length && (!canFind(valid_names, which_half))) {
            throw new Exception("Invalid which_half name: " ~ which_half);
        }
        this.which_half = which_half;
    }

    this(ref const(SpherePatch) other)
    {
        radius = other.radius;
        centre = Vector3(other.centre);
        face_name = other.face_name;
        which_half = other.which_half;
    }

    override SpherePatch dup() const
    {
        return new SpherePatch(radius, centre, face_name, which_half);
    }

    override Vector3 opCall(double r, double s) const
    {
        double x_cube, y_cube, z_cube;
        // First, map the full face of the cube.
        switch (face_name) {
        case "east": x_cube = 1.0; y_cube = -1.0+2.0*r; z_cube = -1.0+2.0*s; break;
        case "west": x_cube = -1.0; y_cube = -1.0+2.0*r; z_cube = -1.0+2.0*s; break;
        case "south": x_cube = -1.0+2.0*r; y_cube = -1.0; z_cube = -1.0+2.0*s; break;
        case "north": x_cube = -1.0+2.0*r; y_cube = 1.0; z_cube = -1.0+2.0*s; break;
        case "bottom": x_cube = -1.0+2.0*r; y_cube = -1.0+2.0*s; z_cube = -1.0; break;
        case "top": x_cube = -1.0+2.0*r; y_cube = -1.0+2.0*s; z_cube = 1.0; break;
        default: throw new Exception("Invalid face: " ~ face_name);
        }
        if (which_half.length > 0) {
            // Limit the mapping to half a face.
            // Note that two of the cube coordinates set above will still be correct.
            // Change just one.
            switch (face_name) {
            case "east":
            case "west":
                switch (which_half) {
                case "south": y_cube = -1.0+r; break;
                case "north": y_cube = r; break;
                case "bottom": z_cube = -1.0+s; break;
                case "top": z_cube = s; break;
                default: throw new Exception("Invalid half-face: " ~ which_half);
                }
                break;
            case "south":
            case "north":
                switch (which_half) {
                case "west": x_cube = -1.0+r; break;
                case "east": x_cube = r; break;
                case "bottom": z_cube = -1.0+s; break;
                case "top": z_cube = s; break;
                default: throw new Exception("Invalid half-face: " ~ which_half);
                }
                break;
            case "bottom":
            case "top":
                switch (which_half) {
                case "west": x_cube = -1.0+r; break;
                case "east": x_cube = r; break;
                case "south": y_cube = -1.0+s; break;
                case "north": y_cube = s; break;
                default: throw new Exception("Invalid half-face: " ~ which_half);
                }
                break;
            default:
                throw new Exception("Invalid face: " ~ face_name);
            }
        }
        double[3] xyz = sphere_mapping(x_cube, y_cube, z_cube);
        return Vector3(xyz[0]*radius+centre.x, xyz[1]*radius+centre.y, xyz[2]*radius+centre.z);
    }

    override string toString() const
    {
        return "SpherePatch(radius=" ~ to!string(radius) ~
            ", centre=" ~ to!string(centre) ~
            ", face=" ~ face_name ~
            ", which_half=" ~ which_half ~
            ")";
    }

private:
    // Map -1,1 cube to a sphere of radius 1.
    double[3] sphere_mapping(double x, double y, double z) const
    {
        double x_dash = x * sqrt(1.0 - 0.5*z*z - 0.5*y*y + y*y*z*z/3.0);
        double y_dash = y * sqrt(1.0 - 0.5*z*z - 0.5*x*x + x*x*z*z/3.0);
        double z_dash = z * sqrt(1.0 - 0.5*y*y - 0.5*x*x + x*x*y*y/3.0);
        return [x_dash, y_dash, z_dash];
    }
    // Source: https://math.stackexchange.com/questions/118760/
    // can-someone-please-explain-the-cube-to-sphere-mapping-formula-to-me
    //
    // For cube with edge length 2a, resulting in sphere with radius a.
    // x' = x * a * math.sqrt( a*a - z*z/2 - y*y/2 + y*y*z*z/(3*a*a))
    // y' = y * a * math.sqrt( a*a - z*z/2 - x*x/2 + x*x*z*z/(3*a*a))
    // z' = z * a * math.sqrt( a*a - y*y/2 - x*x/2 + x*x*y*y/(3*a*a))
} // end class SpherePatch


version(spherepatch_test) {
    import util.msg_service;
    int main() {
        number R = 1.0;
        number zero = 0.0;
        auto east_patch = new SpherePatch(R, Vector3(0.0,0.0,0.0), "east");
        auto mid_e = east_patch(0.5, 0.5);
        auto p5_east = east_patch(0.0, 1.0);
        auto p6_east = east_patch(1.0, 1.0);
        // import std.stdio;
        // writeln("mid_e=", mid_e);
        assert(approxEqualVectors(mid_e, Vector3(R, zero, zero)), failedUnitTest());
        auto top_patch = new SpherePatch(R, Vector3(0.0,0.0,0.0), "top");
        auto mid_t = top_patch(0.5, 0.5);
        auto p5_top = top_patch(1.0, 0.0);
        auto p6_top = top_patch(1.0, 1.0);
        // writeln("mid_t=", mid_t);
        assert(approxEqualVectors(mid_t, Vector3(zero, zero, R)), failedUnitTest());
        assert(approxEqualVectors(p5_top, p5_east), failedUnitTest());
        assert(approxEqualVectors(p6_top, p6_east), failedUnitTest());
        assert(isClose(abs(p6_top), R), failedUnitTest());
        // Now, try a half-face patch.
        auto south_patch = new SpherePatch(R, Vector3(0.0,0.0,0.0), "south", "top");
        auto mid_s = south_patch(0.5, 0.0);
        auto p5_south = south_patch(1.0, 1.0);
        assert(approxEqualVectors(mid_s, Vector3(zero, -R, zero)), failedUnitTest());
        assert(approxEqualVectors(p5_south, p5_east), failedUnitTest());
        return 0;
    }
}

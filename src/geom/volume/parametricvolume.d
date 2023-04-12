/** parametricvolume.d
 * Geometry-building elements for our 3D world -- three-parameter volumes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2015-04-07 first code
 */

module geom.volume.parametricvolume;

import std.stdio;

import geom.elements;

// Nomenclature for the parametric distances, bounding surfaces, paths and corners.
//
// t=1 at top surface
//         north
//    p011-------p111 s=1
//      |         |
// west |   Top   | east
//      |         |
//    p001-------p101 s=0
//         south
//     r=0       r=1
//
//
// t=0 at Bottom surface
//         north
//    p010-------p110 s=1
//      |         |
// west |  Bottom | east
//      |         |
//    p000-------p100 s=0
//         south
//     r=0       r=1
//
// Faces:
// North = face[0]; East = face[1]; South = face[2]; West = face[3];
// Top = face[4]; Bottom = face[5]
//
// Corners:
// Bottom surface: p000 == p[0]; p100 == p[1]; p110 == p[2]; p010 == p[3]
// Top surface   : p001 == p[4]; p101 == p[5]; p111 == p[6]; p011 == p[7]
//
// Edges:
// edge[0] p[0] --> p[1] around Bottom surface
//     [1] p[1] --> p[2]
//     [2] p[3] --> p[2]
//     [3] p[0] --> p[3]
//
//     [4] p[4] --> p[5] around Top surface
//     [5] p[5] --> p[6]
//     [6] p[7] --> p[6]
//     [7] p[4] --> p[7]
//
//     [8] p[0] --> p[4] connecting Bottom to Top
//     [9] p[1] --> p[5]
//    [10] p[2] --> p[6]
//    [11] p[3] --> p[7]
//
// We'll try to use this notation consistently in the classes of this package.

class ParametricVolume {
public:
    abstract Vector3 opCall(double r, double s, double t) const;
    abstract ParametricVolume dup() const;
    abstract override string toString() const;
    Vector3 dVdr(double r, double s, double t) const 
    {
        // Obtain the derivative approximately, via a finite-difference.
        double dr = 0.001;
        Vector3 V0 = this.opCall(r,s,t);
        Vector3 derivative;
        if (r+dr > 1.0) {
            // r is close to the r=1.0 boundary, use a one-sided difference.
            Vector3 Vminus1 = this.opCall(r-dr,s,t);
            derivative = (V0 - Vminus1)/dr;
        } else if (r-dr < 0.0) {
            // r is close to the r=0.0 boundary, use a one-sided difference.
            Vector3 Vplus1 = this.opCall(r+dr,s,t);
            derivative = (Vplus1 - V0)/dr;
        } else {
            // Not near a boundary, use central-difference.
            Vector3 Vminus1 = this.opCall(r-dr,s,t);
            Vector3 Vplus1 = this.opCall(r+dr,s,t);
            derivative = (Vplus1 - Vminus1) / (2.0*dr);
        }
        return derivative;
    }
    Vector3 dVds(double r, double s, double t) const 
    {
        // Obtain the derivative approximately, via a finite-difference.
        double ds = 0.001;
        Vector3 V0 = this.opCall(r,s,t);
        Vector3 derivative;
        if (s+ds > 1.0) {
            // s is close to the s=1.0 boundary, use a one-sided difference.
            Vector3 Vminus1 = this.opCall(r,s-ds,t);
            derivative = (V0 - Vminus1)/ds;
        } else if (s-ds < 0.0) {
            // s is close to the s=0.0 boundary, use a one-sided difference.
            Vector3 Vplus1 = this.opCall(r,s+ds,t);
            derivative = (Vplus1 - V0)/ds;
        } else {
            // Not near a boundary, use central-difference.
            Vector3 Vminus1 = this.opCall(r,s-ds,t);
            Vector3 Vplus1 = this.opCall(r,s+ds,t);
            derivative = (Vplus1 - Vminus1) / (2.0*ds);
        }
        return derivative;
    }
    Vector3 dVdt(double r, double s, double t) const 
    {
        // Obtain the derivative approximately, via a finite-difference.
        double dt = 0.001;
        Vector3 V0 = this.opCall(r,s,t);
        Vector3 derivative;
        if (t+dt > 1.0) {
            // t is close to the t=1.0 boundary, use a one-sided difference.
            Vector3 Vminus1 = this.opCall(r,s,t-dt);
            derivative = (V0 - Vminus1)/dt;
        } else if (t-dt < 0.0) {
            // t is close to the t=0.0 boundary, use a one-sided difference.
            Vector3 Vplus1 = this.opCall(r,s,t+dt);
            derivative = (Vplus1 - V0)/dt;
        } else {
            // Not near a boundary, use central-difference.
            Vector3 Vminus1 = this.opCall(r,s,t-dt);
            Vector3 Vplus1 = this.opCall(r,s,t+dt);
            derivative = (Vplus1 - Vminus1) / (2.0*dt);
        }
        return derivative;
    }
    Vector3 d2Vdr2(double r, double s, double t) const
    {
        // Obtain the derivative approximately, via a finite-difference.
        double dr = 0.001;
        Vector3 V0 = this.opCall(r,s,t);
        Vector3 derivative;
        if (r+dr > 1.0) {
            // r is close to the r=1.0 boundary, use a one-sided difference.
            Vector3 Vminus1 = this.opCall(r-dr,s,t);
            Vector3 Vminus2 = this.opCall(r-2*dr,s,t);
            derivative = (V0 - 2*Vminus1 + Vminus2) / (dr*dr);
        } else if (r-dr < 0.0) {
            // r is close to the r=0 boundary, use a one-sided difference.
            Vector3 Vplus1 = this.opCall(r+dr,s,t);
            Vector3 Vplus2 = this.opCall(r+2*dr,s,t);
            derivative = (Vplus2 - 2*Vplus1 + V0) / (dr*dr);
        } else {
            // Not near a boundary, use central-difference.
            Vector3 Vminus1 = this.opCall(r-dr,s,t);
            Vector3 Vplus1 = this.opCall(r+dr,s,t);
            derivative = (Vplus1 - 2*V0 + Vminus1) / (dr*dr);
        }
        return derivative;
    }
    Vector3 d2Vds2(double r, double s, double t) const
    {
        // Obtain the derivative approximately, via a finite-difference.
        double ds = 0.001;
        Vector3 V0 = this.opCall(r,s,t);
        Vector3 derivative;
        if (s+ds > 1.0) {
            // s is close to the s=1.0 boundary, use a one-sided difference.
            Vector3 Vminus1 = this.opCall(r,s-ds,t);
            Vector3 Vminus2 = this.opCall(r,s-2*ds,t);
            derivative = (V0 - 2*Vminus1 + Vminus2) / (ds*ds);
        } else if (s-ds < 0.0) {
            // s is close to the s=0 boundary, use a one-sided difference.
            Vector3 Vplus1 = this.opCall(r,s+ds,t);
            Vector3 Vplus2 = this.opCall(r,s+2*ds,t);
            derivative = (Vplus2 - 2*Vplus1 + V0) / (ds*ds);
        } else {
            // Not near a boundary, use central-difference.
            Vector3 Vminus1 = this.opCall(r,s-ds,t);
            Vector3 Vplus1 = this.opCall(r,s+ds,t);
            derivative = (Vplus1 - 2*V0 + Vminus1) / (ds*ds);
        }
        return derivative;
    }
    Vector3 d2Vdt2(double r, double s, double t) const
    {
        // Obtain the derivative approximately, via a finite-difference.
        double dt = 0.001;
        Vector3 V0 = this.opCall(r,s,t);
        Vector3 derivative;
        if (t+dt > 1.0) {
            // t is close to the t=1.0 boundary, use a one-sided difference.
            Vector3 Vminus1 = this.opCall(r,s,t-dt);
            Vector3 Vminus2 = this.opCall(r,s,t-2*dt);
            derivative = (V0 - 2*Vminus1 + Vminus2) / (dt*dt);
        } else if (t-dt < 0.0) {
            // t is close to the t=0 boundary, use a one-sided difference.
            Vector3 Vplus1 = this.opCall(r,s,t+dt);
            Vector3 Vplus2 = this.opCall(r,s,t+2*dt);
            derivative = (Vplus2 - 2*Vplus1 + V0) / (dt*dt);
        } else {
            // Not near a boundary, use central-difference.
            Vector3 Vminus1 = this.opCall(r,s,t-dt);
            Vector3 Vplus1 = this.opCall(r,s,t+dt);
            derivative = (Vplus1 - 2*V0 + Vminus1) / (dt*dt);
        }
        return derivative;
    }
    Vector3 d2Vdrds(double r, double s, double t) const {
        // Obtain the derivative approximately, via a finite-difference.
        double dr = 0.001;
        Vector3 dVds0 = this.dVds(r,s,t);
        Vector3 derivative;
        if (r+dr > 1.0) {
            // r is close to the r=1.0 boundary, use a one-sided difference.
            Vector3 dVdsminus1 = this.dVds(r-dr,s,t);
            derivative = (dVds0 - dVdsminus1)/dr;
        } else if (r-dr < 0.0) {
            // r is close to the r=0.0 boundary, use a one-sided difference.
            Vector3 dVdsplus1 = this.dVds(r+dr,s,t);
            derivative = (dVdsplus1 - dVds0)/dr;
        } else {
            // Not near a boundary, use central-difference.
            Vector3 dVdsminus1 = this.dVds(r-dr,s,t);
            Vector3 dVdsplus1 = this.dVds(r+dr,s,t);
            derivative = (dVdsplus1 - dVdsminus1) / (2.0*dr);
        }
        return derivative;
    }
    Vector3 d2Vdrdt(double r, double s, double t) const {
        // Obtain the derivative approximately, via a finite-difference.
        double dr = 0.001;
        Vector3 dVdt0 = this.dVdt(r,s,t);
        Vector3 derivative;
        if (r+dr > 1.0) {
            // r is close to the r=1.0 boundary, use a one-sided difference.
            Vector3 dVdtminus1 = this.dVdt(r-dr,s,t);
            derivative = (dVdt0 - dVdtminus1)/dr;
        } else if (r-dr < 0.0) {
            // r is close to the r=0.0 boundary, use a one-sided difference.
            Vector3 dVdtplus1 = this.dVdt(r+dr,s,t);
            derivative = (dVdtplus1 - dVdt0)/dr;
        } else {
            // Not near a boundary, use central-difference.
            Vector3 dVdtminus1 = this.dVdt(r-dr,s,t);
            Vector3 dVdtplus1 = this.dVdt(r+dr,s,t);
            derivative = (dVdtplus1 - dVdtminus1) / (2.0*dr);
        }
        return derivative;
    }
    Vector3 d2Vdsdt(double r, double s, double t) const {
        // Obtain the derivative approximately, via a finite-difference.
        double ds = 0.001;
        Vector3 dVdt0 = this.dVdt(r,s,t);
        Vector3 derivative;
        if (s+ds > 1.0) {
            // s is close to the s=1.0 boundary, use a one-sided difference.
            Vector3 dVdtminus1 = this.dVdt(r,s-ds,t);
            derivative = (dVdt0 - dVdtminus1)/ds;
        } else if (s-ds < 0.0) {
            // s is close to the s=0.0 boundary, use a one-sided difference.
            Vector3 dVdtplus1 = this.dVdt(r,s+ds,t);
            derivative = (dVdtplus1 - dVdt0)/ds;
        } else {
            // Not near a boundary, use central-difference.
            Vector3 dVdtminus1 = this.dVdt(r,s-ds,t);
            Vector3 dVdtplus1 = this.dVdt(r,s+ds,t);
            derivative = (dVdtplus1 - dVdtminus1) / (2.0*ds);
        }
        return derivative;
    }
        
} // end class ParametricVolume

void writeVolumeAsVtkXml(ParametricVolume vol, string fileName, int nrPts, int nsPts, int ntPts)
{
    double dr = 1.0/(nrPts-1);
    double ds = 1.0/(nsPts-1);
    double dt = 1.0/(ntPts-1);
    auto f = File(fileName, "w");
    f.writeln("<VTKFile type=\"StructuredGrid\" version=\"1.0\" header_type=\"UInt64\">");
    f.writefln("  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">", 0, nrPts-1, 0, nsPts-1, 0, ntPts-1);
    f.writefln("    <Piece Extent=\"%d %d %d %d %d %d\">",  0, nrPts-1, 0, nsPts-1, 0, ntPts-1);
    f.writeln("      <Points>");
    f.writeln("        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">");
    foreach (k; 0 .. ntPts) {
        foreach (j; 0 .. nsPts) {
            foreach (i; 0 .. nrPts) {
                auto r = i*dr;
                auto s = j*ds;
                auto t = k*dt;
                auto p = vol(r, s, t);
                version (complex_numbers) {
                    f.writefln("       %20.16e %20.16e %20.16e", p.x.re, p.y.re, p.z.re);
                } else {
                    f.writefln("       %20.16e %20.16e %20.16e", p.x, p.y, p.z);
                }
            }
        }
    }
    f.writeln("        </DataArray>");
    f.writeln("      </Points>");
    f.writeln("    </Piece>");
    f.writeln("  </StructuredGrid>");
    f.writeln("</VTKFile>");
    f.close();
}

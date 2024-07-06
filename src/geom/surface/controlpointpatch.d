/**
 * An implementation of Control-Point-Form algebraic grid generation
 * for structured grids.
 *
 * Reference:
 * Eiseman (1988)
 * A Control Point Form of Algebraic Grid Generation
 * Internationl Journal of Fluid Mechanics, vol 8, pp 1165--1181
 *
 * Author: Carrie Xie and Rowan G.
 * Date: 2021-12-08
 *
 * Note: Carrie Xie worked on a Python implementation.
 *       Rowan has done a translation and rework into this D implementation.
 */

module geom.surface.controlpointpatch;

import std.conv;
import std.stdio;
import std.format : format;

import geom.elements;
import geom.gpath;
import geom.surface.parametricsurface;
import geom.surface.coonspatch;
import geom.surface.channelpatch;
import geom.surface.aopatch;

class ControlPointPatch : ParametricSurface {
private:
    // Control points are defined in 2D array C[i][j] with
    // 0 <= i <= N-1, and 0 <= j <= M-1
    // Note: we will index from 0 but Eiseman uses 1-based indices
    Vector3[][] mC;
    int mN;
    int mM;
    Path mNorth, mEast, mSouth, mWest;

public:
    this(in Path south, in Path north, in Path west, in Path east, in Vector3[][] C)
    {
        mN = to!int(C.length);
        mM = to!int(C[0].length);

	if (north !is null) mNorth = north.dup();
	if (east !is null) mEast = east.dup();
	if (south !is null) mSouth = south.dup();
	if (west !is null) mWest = west.dup();

	mC.length = mN;
	foreach (i; 0 .. mN) mC[i] = C[i].dup();
    }

    this(in Path south, in Path north, in Path west, in Path east,
         int ncpi, int ncpj, string guidePatch)
    {
        // Let's construct, then alter
        Vector3[][] C;
        C.length = ncpi;
        foreach (ref CC; C) CC.length = ncpj;

	    this(south, north, west, east, C);

        double dr = 1.0/(ncpi - 2);
        double ds = 1.0/(ncpj - 2);
        double[] rPos; rPos.length = ncpi;
        double[] sPos; sPos.length = ncpj;

        rPos[0] = 0.0;
        foreach (i; 1 .. ncpi-1) rPos[i] = dr/2.0 + (i-1)*dr;
        rPos[$-1] = 1.0;

        sPos[0] = 0.0;
        foreach (i; 1 .. ncpj-1) sPos[i] = ds/2.0 + (i-1)*ds;
        sPos[$-1] = 1.0;

        ParametricSurface gpatch;
        bool usingChannelAsGuide = false;
        bool channelBoundedOnEW = false;
        switch (guidePatch) {
        case "coons", "Coons", "COONS", "tfi", "TFI":
            gpatch = new CoonsPatch(south, north, west, east);
            break;
        case "channel", "Channel", "CHANNEL":
            gpatch = new ChannelPatch(south, north, false, true);
            usingChannelAsGuide = true;
            break;
        case "channel-e2w", "Channel-e2w", "CHANNEL-E2W", "channel_e2w":
            gpatch = new ChannelPatch(west, east, false, true);
            usingChannelAsGuide = true;
            channelBoundedOnEW = true;
            break;
        case "aopatch", "AOPatch", "AOPATCH", "area-orthogonal", "area_orthogonal":
            gpatch = new AOPatch(south, north, west, east);
            break;
        default:
            string errMsg = "Error in constructor to ControlPointPatch.\n";
            errMsg ~= format("The guide patch type '%s' is not supported.\n", guidePatch);
            errMsg ~= "Check the spelling or what you are requesting.";
            errMsg ~= "Supported guide patches are: 'coons', 'channel', 'channel-e2w' and 'aopatch'";
            throw new Error(errMsg);
        }

        if (!channelBoundedOnEW) {
            foreach (i; 0 .. ncpi) {
                foreach (j; 0 .. ncpj) {
                    mC[i][j] = gpatch(rPos[i], sPos[j]);
                }
            }
            if (usingChannelAsGuide) {
                // correct the boundaries: east and west
                foreach (j; 0 .. ncpj) {
                    mC[0][j] = west(sPos[j]);
                    mC[$-1][j] = east(sPos[j]);
                }
            }
        }
        else {
            // This case must be with a channel patch bounded on east-west
            foreach (i; 0 .. ncpi) {
                foreach (j; 0 .. ncpj) {
                    mC[i][j] = gpatch(sPos[j], rPos[i]);
                }
            }
            // Correct the boundaries: south and north
            foreach (i; 0 .. ncpi) {
                mC[i][0] = south(rPos[i]);
                mC[i][$-1] = north(rPos[i]);
            }
        }
    }

    this(ref const(ControlPointPatch) other)
    {
	mN = other.mN;
	mM = other.mM;
	if (other.mNorth) mNorth = other.mNorth.dup();
	if (other.mEast) mEast = other.mEast.dup();
	if (other.mSouth) mSouth = other.mSouth.dup();
	if (other.mWest) mWest = other.mWest.dup();
	mC.length = other.mC.length;
	foreach (i; 0 .. mC.length) mC[i] = other.mC[i].dup();
    }

    override ControlPointPatch dup() const
    {
	return new ControlPointPatch(this.mSouth, this.mNorth, this.mWest, this.mEast, this.mC);
    }

    /**
     * From Eisemann's paper, we can selectively apply boundary conformity or not.
     *
     * Eisemann (1998) p. 1172
     * "Variarions of eq. 18 now arise naturally: boundary conformity can be applied selectively.
     * By dropping any conforming term, the associated boundary can be manipulated into the desired shape."
     *
     * In implementation, we'll take the approach of first computing the non-conforming position,
     * and then selectively applying conformity if the edge is specified.
     */
    override Vector3 opCall(double r, double s) const
    {
	// Convert r, s to r_hat, s_hat
	// 0 <= r, s <= 1
	double r_hat = r*(mN-2) + 1.0;
	double s_hat = s*(mM-2) + 1.0;

	Vector3 p = T(r_hat, s_hat);
	if (mWest) p += (1.0 - G(r_hat, 1))*(mWest(s) - F(s_hat, 1));
	if (mEast) p += G(r_hat, mN-1)*(mEast(s) - F(s_hat, mN));
	if (mSouth) p += (1.0 - H(s_hat, 1))*(mSouth(r) - E(r_hat, 1));
	if (mNorth) p += H(s_hat, mM-1)*(mNorth(r) - E(r_hat, mM));
	return p;
    }

    override string toString() const
    {
        return "ControlPointPatch()";
    }

    Vector3 getCtrlPt(int i, int j)
    {
	return mC[i][j];
    }

    void setCtrlPt(int i, int j, Vector3 p)
    {
	mC[i][j] = p;
    }

    void writeCtrlPtsAsVtkXml(string baseFileName)
    {
	auto n = mC.length;
        auto m = mC[0].length;

	// We'll take care of the extension because the VTK-XML files
	// are fussy about extensions when using the XML
	string ext = ".vts";
        auto f = File(baseFileName ~ ext, "w");
        f.writeln("<VTKFile type=\"StructuredGrid\" version=\"1.0\" header_type=\"UInt64\">");
        f.writefln("  <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">", 0, n-1, 0, m-1);
        f.writefln("    <Piece Extent=\"%d %d %d %d 0 0\">",  0, n-1, 0, m-1);
        f.writeln("      <Points>");
        f.writeln("        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">");
        foreach (j; 0 .. m) {
            foreach (i; 0 .. n) {
                auto p = mC[i][j];
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
private:
    double Omega(double r) const
    {
	// Implements Eq (4) in Eiseman (1988)
	if (r < -1.0) return 0.0;
	if (-1.0 <= r && r < 0.0) return 0.5*(r + 1.0)*(r + 1.0);
	if (0 <= r && r <= 1.0) return 1.0 - 0.5*(r - 1.0)*(r - 1.0);
	// all other cases, 1 < r
	return 1.0;
    }

    double G(double r, int alpha) const
    {
	// Implements Eq (6) in Eiseman (1988)
	// Deal with special end cases
	if (alpha == 1) return 2.0*Omega(r-1.0) - 1.0;
	if (alpha == mN-1) return 2.0*Omega(r - mN + 1.0);
	// For other cases
	return Omega(r - alpha);
    }

    double H(double s, int beta) const
    {
        if (beta == 1) return 2.0*Omega(s-1.0) - 1.0;
	if (beta == mM-1) return 2.0*Omega(s - mM + 1.0);
	// For other cases
	return Omega(s - beta);
    }


    Vector3 E(double r, int j) const
    {
	// REVISIT A POTENTIAL OPTIMISATION
	// It looks like we only ever use the displacement Vectors (C_alpha,j - C_alpha-1,j)
	// so we could pre-compute and store these.
	// (Plus we'd need starting points.)

	// Implements Eq (7) in Eiseman (1988)
	// The alpha index is annoying here.
	// We need to work with 1-based in the G function call.
	// We need to work with 0-based in the access to mC.
	// This is why C_{alpha+1,j} in Eiseman's equation becomes: mC[alpha][j]
	Vector3 p = mC[0][j-1];
	foreach (alpha; 1 .. mN) {
	    p += G(r, alpha)*(mC[alpha][j-1] - mC[alpha-1][j-1]);
	}
	return p;
    }

    Vector3 F(double s, int i) const
    {
	// Implements Eq (8) in Eismeann (1988)
	// As above for E function, now the beta index is annoying.
	// We need to work with 1-based in the H function call.
	// We need to work with 0-based in the access to mC.
	// This is why C_{i,beta+1} in Eiseman's equation becomes: mC[i][beta]
	Vector3 p = mC[i-1][0];
	foreach (beta; 1 .. mM) {
	    p += H(s, beta)*(mC[i-1][beta] - mC[i-1][beta-1]);
	}
	return p;
    }

    Vector3 T(double r, double s) const
    {
	// Implements Eq (9) in Eiseman (1988)
	Vector3 p = E(r, 1);
	foreach (beta; 1 .. mM) {
	    p += H(s, beta)*(E(r, beta+1) - E(r, beta));
	}
	return p;
    }

}

version(controlpointpatch_test) {

    import util.msg_service;
    int main() {
        int N = 4;
        int M = 5;
        double L = 1.0;

        Vector3[][] C;
        C.length = N;
        foreach (i; 0 .. N) C[i].length = M;

        C[0][0] = Vector3(0.0, 0.0);
        C[1][0] = Vector3(0.25*L, 0.0);
        C[2][0] = Vector3(0.75*L, 0.0);
        C[3][0] = Vector3(L, 0.0);

        auto dy = Vector3(0.0, L/6.0);
        foreach (i; 0 .. N) C[i][1] = C[i][0] + dy;

        dy = Vector3(0.0, 3.0*L/6.0);
        foreach (i; 0 .. N) C[i][2] = C[i][0] + dy;

        dy = Vector3(0.0, 5.0*L/6.0);
        foreach (i; 0 .. N) C[i][3] = C[i][0] + dy;

        dy = Vector3(0.0, L);
        foreach (i; 0 .. N) C[i][4] = C[i][0] + dy;

        Line north = new Line(C[0][4], C[3][4]);
        Line east = new Line(C[3][0], C[3][4]);
        Line south = new Line(C[0][0], C[3][0]);
        Line west = new Line(C[0][0], C[0][4]);

        ControlPointPatch ctrlPtPatch = new ControlPointPatch(south, north, west, east, C);

        auto p = ctrlPtPatch(0.5, 0.5);
        assert(approxEqualVectors(p, Vector3(0.5, 0.5)), failedUnitTest());
        p = ctrlPtPatch(1.0, 1.0);
        assert(approxEqualVectors(p, Vector3(1.0, 1.0)), failedUnitTest());

        return 0;
    }
}

/** nurbsvolume.d
 * Author: Reece O. 
 * Date: 2021-03-09
 */

module geom.volume.nurbsvolume;

import std.array;
import std.conv;
import std.format;
import std.file;
import std.stdio;
import std.string;

import geom.elements;
import geom.volume.parametricvolume;
import geom.misc.nurbs_utils;

import nm.number;
import ntypes.complex;

class NURBSVolume : ParametricVolume {
public:
    number[4][][][] Pw; // weighted control point array
    double[] U, V, W;   // knot vectors
    int p, q, r;        // spline degrees
    int nPu, nPv, nPw;  // dimensions of control net
    int nU, nV, nW;     // lengths of knot vector arrays

    this(const number[4][][][] Pw, const double[] U, const int p, const double[] V, 
            const int q, const double[] W, const int r) 
    {
        initialise(Pw, U, p, V, q, W, r);
    }
    this(string fileName, string fmt, int[3] degrees, double scale=1.0){
        readFromCntrlNetFile(fileName, fmt, degrees, scale);
    }
    this(string fileName, double scale=1.0){
        readFromFile(fileName, scale); 
    }
    this(ref const NURBSVolume other)
    {
        this.Pw = duplicatePw(other.Pw); 
        this.p = other.p; 
        this.q = other.q; 
        this.r = other.r;
        this.U = other.U.dup; 
        this.V = other.V.dup; 
        this.W = other.W.dup; 
        this.nPu = other.nPu;
        this.nPv = other.nPv;
        this.nPw = other.nPw;
        this.nU = other.nU;
        this.nV = other.nV;
        this.nW = other.nW;
        this.m_a = other.m_a;
        this.m_b = other.m_b;
        this.m_c = other.m_c;
        mNu.length = other.p + 1;
        mNws_u = NURBSWorkspace(other.p);
        mNv.length = other.q + 1;
        mNws_v = NURBSWorkspace(other.q);
        mNw.length = other.r + 1;
        mNws_w = NURBSWorkspace(other.r);
        temp1.length = q+1;
        foreach (i; 0 .. temp1.length) temp1[i].length = r+1;
	    temp2.length = r+1;
    }
    override NURBSVolume dup() const
    {
        return new NURBSVolume(Pw, U, p, V, q, W, r);
    }
    override Vector3 opCall(double u, double v, double w) const
    {
        return eval(u, v, w);
    }
    override string toString() const
    {
        string repString = "NURBSVolume(Pw=";
        repString ~= to!string(Pw);
        repString ~= ", U=" ~ to!string(U) ~ ", p=" ~ to!string(p);
        repString ~= ", V=" ~ to!string(V) ~ ", q=" ~ to!string(q);
        repString ~= ", W=" ~ to!string(W) ~ ", r=" ~ to!string(r) ~ ")";
        return repString;
    }
    
    void initialise(const number[4][][][] Pw, const double[] U, int p, const double[] V, 
            int q, const double[] W, int r)
    {
        // Test Pw is minimum viable
        nPu = to!int(Pw.length);
        nPv = to!int(Pw[0].length);
        nPw = to!int(Pw[0][0].length);
	if (nPu < 2 || nPv < 2) {
	    string errMsg = "NURBSVolume() A NURBS volume requires at least two ";
            errMsg ~= "control points in each direction.\n";
	    errMsg ~= format("Supplied number of control points in direction u: %s\n", nPu);
	    errMsg ~= format("Supplied number of control points in direction v: %s", nPv);
	    throw new Error(text(errMsg));
	}
	
        // Test that p is valid
        m_a = to!int(nPu-1);
	if (p > m_a) {
	    string errMsg = "NURBSVolume() Degree of NURBS volume in direction u is not ";
            errMsg ~= "compatible with given number of control points.\n";
	    errMsg ~= format("Supplied degree in u direction: %s\n", p);
	    errMsg ~= format("Maximum allowed degree in u direction: %s", Pw.length-1);
	    throw new Error(text(errMsg));
	}
	
        // Test that q is valid
        m_b = to!int(nPv-1);
	if (q > m_b) {
	    string errMsg = "NURBSVolume() Degree of NURBS volume in direction v is not ";
            errMsg ~= "compatible with given number of control points.\n";
	    errMsg ~= format("Supplied degree in v direction: %s\n", q);
	    errMsg ~= format("Maximum allowed degree in v direction: %s", Pw[0].length-1);
	    throw new Error(text(errMsg));
	}
	
        // Test that r is valid
        m_c = to!int(nPw-1); 
	if (r > m_c) {
	    string errMsg = "NURBSVolume() Degree of NURBS volume in direction w is not ";
            errMsg ~= "compatible with given number of control points.\n";
	    errMsg ~= format("Supplied degree in w direction: %s\n", r);
	    errMsg ~= format("Maximum allowed degree in w direction: %s", Pw[0][0].length-1);
	    throw new Error(text(errMsg));
	}
	
        // Test U knot vector is valid
        nU = to!int(U.length); 
	if (nU != m_a+p+2) {
            string errMsg = "NURBSVolume() Knot vector in direction u is not compatible ";
            errMsg ~= "with given number control points and volume degree.\n";
            errMsg ~= format("Supplied number of knots in direction u: %s\n", U.length);
            errMsg ~= format("Required number of knots in direction u: %s", m_a+p+2);
            throw new Error(text(errMsg));
        }
        
        // Test V knot vector is valid
	nV = to!int(V.length);
        if (V.length != m_b+q+2) {
            string errMsg = "NURBSVolume() Knot vector in direction v is not compatible ";
            errMsg ~= "with given number control points and volume degree.\n";
            errMsg ~= format("Supplied number of knots in direction v: %s\n", V.length);
            errMsg ~= format("Required number of knots in direction v: %s", m_b+q+2);
            throw new Error(text(errMsg));
        }
        
        // Test W knot vector is valid
	nW = to!int(W.length);
        if (W.length != m_c+r+2) {
            string errMsg = "NURBSVolume() Knot vector in direction w is not compatible ";
            errMsg ~= "with given number control points and volume degree.\n";
            errMsg ~= format("Supplied number of knots in direction w: %s\n", W.length);
            errMsg ~= format("Required number of knots in direction w: %s", m_c+r+2);
            throw new Error(text(errMsg));
        }
        
        this.Pw = duplicatePw(Pw);
        this.p = p; 
        this.q = q; 
        this.r = r;
        this.U = U.dup; 
        this.V = V.dup; 
        this.W = W.dup; 
        this.nPu = nPu;
        this.nPv = nPv;
        this.nPw = nPw;
        this.nU = nU;
        this.nV = nV;
        this.nW = nW;
        this.m_a = m_a;
        this.m_b = m_b;
        this.m_c = m_c;
        mNu.length = p+1;
        mNws_u = NURBSWorkspace(p);
        mNv.length = q+1;
        mNws_v = NURBSWorkspace(q);
        mNw.length = r+1;
        mNws_w = NURBSWorkspace(r);
        temp1.length = q+1;
        foreach (i; 0 .. temp1.length) temp1[i].length = r+1;
	    temp2.length = r+1;
    }

    void readFromCntrlNetFile(string fileName, string fmt, int[3] degrees, double scale) 
    {
        /**
         * Constructs NURBS volume from control net file.
         * 
         * Supported file formats:
         *   "plot3d-sgrid": Plot3D structured grid format (order n,k,j,i)
         *   "sgrid": standard structured grid format (order i,j,k,n)
         */
        
        // import hull file
        auto hullData = readText(fileName).split;
        
        // get number of control points in each direction
        int nPu = to!int(hullData[0]);
        int nPv = to!int(hullData[1]);
        int nPw = to!int(hullData[2]);
        
        // initialise control point array
        number[4][][][] Pw = new number[4][][][](nPu, nPv, nPw);
        int ni, nj, nk, nl;
        int iOffset, jOffset, kOffset, offset;
        int headOffset = 3;
    
        // determine order in which control points are imported
        if (fmt == "plot3d-sgrid") {
            ni = 3; nj = nPw; nk = nPv; nl = nPu;
        } else if (fmt == "sgrid") {
            ni = nPu; nj = nPv; nk = nPw; nl = 3;
        } else {
            string errMsg = format("Invalid NURBS net format: '%s'. ", fmt);
            errMsg ~= "Accepted options are 'plot3d-sgrid' or 'sgrid'.";
            throw new Error(errMsg);
        }

        // import control net
        foreach (i; 0 .. ni) {
            iOffset = i*nj*nk*nl;
            foreach (j; 0 .. nj) {
                jOffset = j*nk*nl;
                foreach (k; 0 .. nk) {
                    kOffset = k*nl;
                    foreach (n; 0 .. nl) {
                        // add header offset so cntrl net dimensions are not included
                        offset = iOffset + jOffset + kOffset + headOffset + n; 
                        if (fmt == "plot3d-sgrid") {
                            Pw[n][k][j][i] = to!number(hullData[offset])*scale;
                        } else if (fmt == "sgrid") {
                            Pw[i][j][k][n] = to!number(hullData[offset])*scale;
                        }
                    }
                }
            }
        }
    
        // set all weights equal to 1.0 for initial hull
        foreach (i; 0 .. nPu) {
            foreach (j; 0 .. nPv) {
                foreach (k; 0 .. nPw) {
                    Pw[i][j][k][3] = 1.0;
                }
            }
        }

        // create uniform knot vectors that are clamped and normalised
        int p = degrees[0];
        int q = degrees[1];
        int r = degrees[2];
        auto U = autoKnotVector(nPu, p);
        auto V = autoKnotVector(nPv, q);
        auto W = autoKnotVector(nPw, r);
        
        // initialise NURBS volume
        initialise(Pw, U, p, V, q, W, r);
    }
    
    void readFromFile(string fileName, double scale=1.0)
    {
        int p, q, r, nU, nV, nW, nPu, nPv, nPw;
        double[] U, V, W;

        string[] tokens;
        string token;
        auto f = File(fileName, "r");

        // degrees
        tokens = f.readln().strip().split();
        p = to!int(tokens[0]);
        q = to!int(tokens[1]);
        r = to!int(tokens[2]);

        // knot vector U
        token = f.readln().strip();
        nU = to!int(token);
        U.length = nU;
        foreach(i; 0 .. nU){
            token = f.readln().strip();
            U[i] = to!double(token);
        }

        // knot vector V
        token = f.readln().strip();
        nV = to!int(token);
        V.length = nV;
        foreach(i; 0 .. nV){
            token = f.readln().strip();
            V[i] = to!double(token);
        }

        // knot vector W
        token = f.readln().strip();
        nW = to!int(token);
        W.length = nW;
        foreach(i; 0 .. nW){
            token = f.readln().strip();
            W[i] = to!double(token);
        }

        // control net dimensions
        tokens = f.readln().strip().split();
        nPu = to!int(tokens[0]);
        nPv = to!int(tokens[1]);
        nPw = to!int(tokens[2]);

        // import control net
        number[4][][][] Pw = new number[4][][][](nPu, nPv, nPw);
        foreach(i; 0 .. nPu) {
            foreach(j; 0 .. nPv) {
                foreach(k; 0 .. nPw) {
                    tokens = f.readln().strip().split();
                    foreach(n; 0 .. 4) { 
                        Pw[i][j][k][n] = to!double(tokens[n]);
                    }
                }
            }   
        }

        // multiply control point coords by respective weight and scale
        foreach(i; 0 .. nPu) {
            foreach(j; 0 .. nPv) {
                foreach(k; 0 .. nPw) {
                    foreach(n; 0 .. 3) {
                        Pw[i][j][k][n] *= scale*Pw[i][j][k][3];
                    }
                }
            }
        }

        // initialise NURBS volume
        initialise(Pw, U, p, V, q, W, r);
    }
    
    void writeToFile(string fileName)
    {
        auto f = File(fileName, "w");
        string line;

        f.writefln("%d %d %d", this.p, this.q, this.r);
        f.writefln("%d", this.nU);
        foreach(i; 0 .. this.nU) {
            f.writefln("%.18e", this.U[i]);
        }
        f.writefln("%d", this.nV);
        foreach(i; 0 .. this.nV) {
            f.writefln("%.18e", this.V[i]);
        }
        f.writefln("%d", this.nW);
        foreach(i; 0 .. this.nW) {
            f.writefln("%.18e", this.W[i]);
        } 
        f.writefln("%d %d %d", this.nPu, this.nPv, this.nPw);
        foreach(i; 0 .. nPu) {
            foreach(j; 0 .. nPv) {
                foreach(k; 0 .. nPw) {
                    line = "";
                    foreach(n; 0 .. 3) { 
                        line ~= format("%.18e ", 
                                this.Pw[i][j][k][n].re/this.Pw[i][j][k][3].re); 
                    }
                    line ~= format("%.18e", this.Pw[i][j][k][3].re);
                    f.writeln(line);
                }
            }
        } 
    }
    
    void writeAsVtkXml(string fileName, int nuPts=20, int nvPts=20, int nwPts=20)
    {
        writeVolumeAsVtkXml(this, fileName, nuPts, nvPts, nwPts);
    }

    void writeCntrlNetAsVtkXml(string fileName)
    {
        // TODO: just convert control net as sgrid, then export
        Vector3 p;
        int nrPts = this.nPu;
        int nsPts = this.nPv;
        int ntPts = this.nPw;
                
        auto f = File(fileName, "w");
        f.writeln("<VTKFile type=\"StructuredGrid\" version=\"1.0\" header_type=\"UInt64\">");
        f.writefln("  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">", 0, nrPts-1, 0, nsPts-1, 0, ntPts-1);
        f.writefln("    <Piece Extent=\"%d %d %d %d %d %d\">",  0, nrPts-1, 0, nsPts-1, 0, ntPts-1);
        f.writeln("      <Points>");
        f.writeln("        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">");
        foreach (k; 0 .. ntPts) {
            foreach (j; 0 .. nsPts) {
                foreach (i; 0 .. nrPts) {
                    p.x = this.Pw[i][j][k][0]/this.Pw[i][j][k][3];
                    p.y = this.Pw[i][j][k][1]/this.Pw[i][j][k][3];
                    p.z = this.Pw[i][j][k][2]/this.Pw[i][j][k][3];
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
    
private:
    int m_a; // a+1 = number of control points in direction u
    int m_b; // b+1 = number of control points in direction v
    int m_c; // c+1 = number of control points in direction w
    
    static number[4] mSw;
    static number[3] mS; 
    static number[4][][] temp1; 
    static number[4][] temp2;
    static double[] mNu;
    static double[] mNv;
    static double[] mNw;
    static NURBSWorkspace mNws_u;
    static NURBSWorkspace mNws_v;
    static NURBSWorkspace mNws_w;
    
    Vector3 eval(double u, double v, double w) const {
        // Returns the Cartesian coordinates of a point on a NURBS volume at a 
        // given set of parameter values
        int uSpan = findSpan(u, m_a, p, U);
        basisFuns(uSpan, u, p, U, mNu, mNws_u);
        int vSpan = findSpan(v, m_b, q, V);
        basisFuns(vSpan, v, q, V, mNv, mNws_v);
        int wSpan = findSpan(w, m_c, r, W);
        basisFuns(wSpan, w, r, W, mNw, mNws_w);
        foreach (m; 0 .. r+1) {
            foreach (l; 0 .. q+1) {
                foreach (n; 0 .. 4) { temp1[l][m][n] = 0.0; }
                foreach (k; 0 .. p+1) {
                    foreach (n; 0 .. 4) {
                        temp1[l][m][n] += mNu[k]*Pw[uSpan-p+k][vSpan-q+l][wSpan-r+m][n];
                    }
                }
            }
        }
        foreach (m; 0 .. r+1) {
            foreach (n; 0 .. 4) { temp2[m][n] = 0.0; }
            foreach (l; 0 .. q+1) {
                foreach (n; 0 .. 4) {
                    temp2[m][n] += mNv[l]*temp1[l][m][n];
                }
            }
        }
        foreach (n; 0 .. 4) { mSw[n] = 0.0; }
        foreach (n; 0 .. 3) { mS[n] = 0.0; }
        foreach (m; 0 .. r+1) {
            foreach (n; 0 .. 4) {
                mSw[n] += mNw[m]*temp2[m][n];
            }
        }
        foreach (i; 0 .. 3) { mS[i] = mSw[i]/mSw[3]; }
        return Vector3(mS);
    }
}

version(nurbsvolume_test) {
    import util.msg_service;
    import std.math;
    int main () {
        // solid cylinder point evaluation test
        double s2 = sqrt(2.0);
        number[4][][][] Pw = [[[
            [0.0, 0.0, -3.0, 1.0], [0.0, 0.0, -3.0, 1.0], [0.0, 0.0, -3.0, 1.0], 
            [0.0, 0.0, -3.0, 1.0], [0.0, 0.0, -3.0, 1.0], [0.0, 0.0, -3.0, 1.0],             
            [0.0, 0.0, -3.0, 1.0],  [0.0, 0.0, -3.0, 1.0], [0.0, 0.0, -3.0, 1.0]],
                
            [[0.0, 1.0, -3.0, 1.0], [1.0, 1.0, -3.0, 1.0/s2], [1.0, 0.0, -3.0, 1.0], 
            [1.0, -1.0, -3.0, 1.0/s2], [0.0, -1.0, -3.0, 1.0], [-1.0, -1.0, -3.0, 1.0/s2], 
            [-1.0, 0.0, -3.0, 1.0], [-1.0, 1.0, -3.0, 1.0/s2], [0.0, 1.0, -3.0, 1.0]]],
                              
            [[[0.0, 0.0, -1.0, 1.0], [0.0, 0.0, -1.0, 1.0], [0.0, 0.0, -1.0, 1.0], 
            [0.0, 0.0, -1.0, 1.0], [0.0, 0.0, -1.0, 1.0], [0.0, 0.0, -1.0, 1.0],             
            [0.0, 0.0, -1.0, 1.0],  [0.0, 0.0, -1.0, 1.0],  [0.0, 0.0, -1.0, 1.0]],

            [[0.0, 1.0, -1.0, 1.0], [1.0, 1.0, -1.0, 1.0/s2], [1.0, 0.0, -1.0, 1.0], 
            [1.0, -1.0, -1.0, 1.0/s2], [0.0, -1.0, -1.0, 1.0], [-1.0, -1.0, -1.0, 1.0/s2], 
            [-1.0, 0.0, -1.0, 1.0], [-1.0, 1.0, -1.0, 1.0/s2], [0.0, 1.0, -1.0, 1.0]]],
                               
            [[[0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0], 
            [0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0],             
            [0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0]],
                
            [[0.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0/s2], [1.0, 0.0, 1.0, 1.0], 
            [1.0, -1.0, 1.0, 1.0/s2], [0.0, -1.0, 1.0, 1.0], [-1.0, -1.0, 1.0, 1.0/s2], 
            [-1.0, 0.0, 1.0, 1.0], [-1.0, 1.0, 1.0, 1.0/s2], [0.0, 1.0, 1.0, 1.0]]],
                               
            [[[0.0, 0.0, 3.0, 1.0], [0.0, 0.0, 3.0, 1.0], [0.0, 0.0, 3.0, 1.0], 
            [0.0, 0.0, 3.0, 1.0], [0.0, 0.0, 3.0, 1.0],  [0.0, 0.0, 3.0, 1.0],             
            [0.0, 0.0, 3.0, 1.0],  [0.0, 0.0, 3.0, 1.0], [0.0, 0.0, 3.0, 1.0]],
                
            [[0.0, 1.0, 3.0, 1.0], [1.0, 1.0, 3.0, 1.0/s2], [1.0, 0.0, 3.0, 1.0], 
            [1.0, -1.0, 3.0, 1.0/s2], [0.0, -1.0, 3.0, 1.0], [-1.0, -1.0, 3.0, 1.0/s2], 
            [-1.0, 0.0, 3.0, 1.0], [-1.0, 1.0, 3.0, 1.0/s2], [0.0, 1.0, 3.0, 1.0]]]];
        
        int p = 3;
        int q = 1;
        int r = 2;
        double[] Z = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        double[] R = [0.0, 0.0, 1.0, 1.0];
        double[] T = [0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0];
        
        auto nvol = new NURBSVolume(Pw, Z, p, R, q, T, r);
        Vector3 P = Vector3(0.5, 0.0, 1.5);
        assert(approxEqualVectors(P, nvol(0.75, 0.5, 0.25)), failedUnitTest());

        return 0;
    }
}

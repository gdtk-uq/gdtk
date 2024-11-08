/** 
 * nurbs_utils.d
 * Author: Reece O. 
 * Date: 2021-03-09
 */

module geom.misc.nurbs_utils;

import std.conv;
import std.format;
import std.math;

import geom.volume.nurbsvolume;

import nm.number;

int findSpan(double u, int n, int p, const double[] U) {
    // Returns the index of the given knot vector whose value is less than the given parameter
    // This is algorithm A2.1 from Piegl and Tiller (1997) - 'The NURBS Book'
    
    // check that given parameter is within given parameter range
    // adjust parameter value if slightly out of range 
    double tol = 1.0e-12;
    if ((u < U[0]) && (U[0]-u < tol)) {
        u = U[0];
    } else if ((u > U[$-1]) && (u-U[$-1] < tol)) {
        u = U[$-1];
    } else if (u < U[0] || u > U[$-1]) {
        string errMsg = "Error in FindSpan.\n";
        errMsg ~= "The parameter value is not compatible with the chosen knot vector.\n";
        errMsg ~= format("Supplied parameter value: %.18e\n", u);
        errMsg ~= format("Valid parameter range: [%.18e, %.18e]", U[0], U[U.length-1]);
        throw new Exception(errMsg);
    }
    
    // computing span index
    if (u == U[n+1]) return n;
    auto low = p;
    auto high = n+1;
    auto mid = (low + high)/2;
    while (u < U[mid] || u >= U[mid+1]) {
        if (u < U[mid]) high = mid;
        else low = mid;
        mid = (low+high)/2;
    }
    return mid;
}

struct NURBSWorkspace {
    double[] left;
    double[] right;
    this(int p) {
        this.left.length = p+1;
        this.right.length = p+1; 
    }
}

void basisFuns(int i, double u, int p, const double[] U, ref double[] N, ref NURBSWorkspace nws) {    
    // Fills in N array with all nonvanishing basis function terms
    // This is algorithm A2.2 from Piegl and Tiller (1997) - 'The NURBS Book'
    N[0] = 1.0;
    foreach (j; 1 .. p+1) {
        nws.left[j] = u-U[i+1-j];
        nws.right[j] = U[i+j]-u;
        double saved = 0.0;
        foreach (r; 0 .. j) {
            auto temp = N[r]/(nws.right[r+1] + nws.left[j-r]);
            N[r] = saved + nws.right[r+1]*temp;
            saved = nws.left[j-r]*temp;
        }
        N[j] = saved;
    }
}

void PwTest(const double[4][] Pw) {
    if (Pw.length < 2) {
        string errMsg = "NURBS() A NURBS path requires at least two control points.\n";
        errMsg ~= format("Supplied number of control points: %s", Pw.length);
        throw new Error(text(errMsg));
    }
}

number[4][][][] duplicatePw(const number[4][][][] Pw) {
    number[4][][][] Pw_copy;
    Pw_copy.length = Pw.length;
    foreach (i, array_i; Pw) {
        Pw_copy[i].length = array_i.length;
        foreach (j, array_j; array_i) {
            Pw_copy[i][j].length = array_j.length;
            foreach (k, array_k; array_j) {
                Pw_copy[i][j][k] = array_k;
            }
        }
    }
    return Pw_copy;
}

double[] autoKnotVector(int N, int p) {
    // creates a knot vector this is clamped, uniform and normalised
    // ensure number of control points and degree are compatible
    int a = N - 1;
    if ((p > a) || (p < 1)) {
        throw new Error("Number of control points is not compatible with degree (1<=p<=N+1).");
    }
    int q = a + p + 1;
    double[] U;
    U.length = q + 1;
    foreach (i; 0 .. q+1) {
        // clamp start of curve
        if ((0 <= i) && (i <= p)) { U[i] = 0.0; }
        // fill internal knots
        if ((p < i) && (i <= q-p-1)) { U[i] = i - p; }
        // clamp end of curve
        if ((q-p-1 < i) && (i <= q)) { U[i] = q - 2*p; }
    }
    // normalise knot vector
    double UMax = U[$-1];
    foreach (i; 0 .. q+1) { U[i] /= UMax; }
    return U;
}

void convNURBSVolDataToVTK(string nurbsDataFile, string volVTKFile="nurbs.vts", string netVTKFile="net.vts")
{
    auto nurbs = new NURBSVolume(nurbsDataFile);
    nurbs.writeAsVtkXml(volVTKFile);
    nurbs.writeCntrlNetAsVtkXml(netVTKFile);
    return;
}

version(nurbs_utils_test) {
    import util.msg_service;
    int main() {
        // example 2.3 from Piegl and Tiller (1997) - 'The NURBS Book'
        // findSpan test
        double u = 2.5;
        double[] U = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0];
        int m = 11;
        int p = 2;
        int n = m - p - 1;
        int i = findSpan(u, n, p, U);
        int iExact = 4;
        assert(isClose(i, iExact), failedUnitTest());
        
        // basisFuns test
        auto nws = NURBSWorkspace(p);
        double[] N;
        N.length = p + 1;
        basisFuns(i, u, p, U, N, nws);
        double[3] Nexact = [1.0/8.0, 3.0/4.0, 1.0/8.0];
        assert((N.length == Nexact.length), failedUnitTest());
        foreach (idx; 0 .. N.length) {
            assert(isClose(N[idx], Nexact[idx]), failedUnitTest());
        }

        // single-span autoKnotVector test
        int p_ss = 3;
        int N_ss = 4;
        auto U_ss = autoKnotVector(N_ss, p_ss);
        double[8] UExact_ss = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        assert((U_ss.length == UExact_ss.length), failedUnitTest());
        foreach(idx; 0 .. U_ss.length) {
            assert((isClose(U_ss[idx], UExact_ss[idx])), failedUnitTest());
        }
        
        // multi-span autoKnotVector test
        int p_ms = 2;
        int N_ms = 5;
        auto U_ms = autoKnotVector(N_ms, p_ms);
        double[8] UExact_ms = [0.0, 0.0, 0.0, 1.0/3.0, 2.0/3.0, 1.0, 1.0, 1.0];
        assert((U_ms.length == UExact_ms.length), failedUnitTest());
        foreach(idx; 0 .. U_ms.length) {
            assert((isClose(U_ms[idx], UExact_ms[idx])), failedUnitTest());
        }
        
        return 0;
    }
}

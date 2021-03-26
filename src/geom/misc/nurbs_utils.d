/** 
 * nurbs_utils.d
 * Author: Reece O. 
 * Date: 2021-03-09
 */

module nurbs_utils;

import std.format;
import std.conv;

int FindSpan(double u, int n, int p, const double[] U) {
    // Returns the index of the given knot vector whose value is less than the given parameter
    // This is algorithm A2.1 from Piegl and Tiller (1997) - 'The NURBS Book'
    
    // check that given parameter is within given parameter range
    if (u < U[0] || u > U[U.length-1]) {
        string errMsg = "Error in FindSpan: the supplied parameter value is not within the supplied parameter range.\n";
        errMsg ~= format("Supplied parameter value: %f\n", u);
        errMsg ~= format("Supplied parameter range: [%(%f, %)]", [U[0], U[U.length-1]]);
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

void BasisFuns(int i, double u, int p, const double[] U, ref double[] N, ref NURBSWorkspace nws) {    
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



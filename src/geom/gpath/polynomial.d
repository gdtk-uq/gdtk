// polynomial.d

module geom.gpath.polynomial;

import std.conv;
import std.math;

import nm.bbla;
import ntypes.complex;
import nm.number;

import geom.elements;
import geom.gpath.path;


class Polynomial : Path {
    // polynomial in the form y = sum(c_i * x^^i)|i=0..n
    // Momar Hughes, 4th-year thesis, October 2015
public:
    Vector3[] P; // array of control points to interpolate
    number[] C; // array of coefficients

    this(in Vector3[] P)
    {
        if (P.length == 0) {
            throw new Error(text("Polynomial() No control points present."));
        } // end if
        if (P.length == 1) {
            throw new Error(text("Polynomial() Only one control point, not enough for a curve."));
        } // end if
        this.P = P.dup();
        evaluate_coefficients();
    }
    this(in number[] C, number x0, number x1)
    {
        if (C.length == 0) {
            throw new Error(text("Polynomial() No coefficients provided."));
        } // end if
        this.C = C.dup();
        this.P.length = 2;
        this.P[0] = evaluate_polynomial(x0);
        this.P[1] = evaluate_polynomial(x1);
    }
    this(in Vector3[] P,in number[] C)
    {
        if (C.length == 0) {
            throw new Error(text("Polynomial() No coefficients provided."));
        } // end if
        this.C = C.dup();
        foreach(point;P){
            if(point.y != evaluate_polynomial(point.x).y){
                throw new Error(text("Polynomial() points and coefficients do not match."));
            } // end if
        } // end foreach
        this.P = P.dup();
    }
    this(ref const(Polynomial) other)
    {
        this.P = other.P.dup();
        this.C = other.C.dup();
    }
    override Polynomial dup() const
    {
        return new Polynomial(P,C);
    } // end dup()
    override Vector3 opCall(double t) const
    {
        // Evaluate P(t)
        number xt = P[0].x + t*(P[$-1].x-P[0].x);
        return evaluate_polynomial(xt);
    } // end opCall()
    override Vector3 dpdt(double t) const
    {
        // Evaluate P(t)
        number xt;
        xt = P[0].x + t*(P[$-1].x-P[0].x);
        return derivative_polynomial(xt);
    }
    override string toString() const
    {
        return "Polynomial(P=" ~ to!string(P) ~ ")";
    }
    override string classString() const
    {
        return "Polynomial";
    }

protected:
    void evaluate_coefficients()
    {
        size_t n = P.length;
        auto A = new Matrix!number(n);
        auto b = new Matrix!number(n,1);
        foreach(i;0 .. n){
            b[i,0] = P[i].y;
            foreach(j;0 .. n){
                A[i,j] = P[i].x^^j;
            } // end foreach
        } // end foreach
        auto x = lsqsolve!number(A,b);
        foreach(i;0 .. n){
            C ~= x[i,0];
        } // end foreach
    } // end evaluate_coefficients ()

    Vector3 evaluate_polynomial(number x) const
    {
        number y=0.0;
        foreach(i;0 .. C.length){
            y += C[i] * x^^i;
        } // end foreach
        return Vector3(x,y);
    } // end evaluate_polynomial ()

    Vector3 derivative_polynomial(number x) const
    {
        number dy=0.0;
        foreach(i;0 .. C.length){
            dy += C[i] * i * x^^(i-1);
        } // end foreach
        return Vector3(to!number(1.0),dy);
    } // end evaluate_polynomial ()

} // end class Polynomial

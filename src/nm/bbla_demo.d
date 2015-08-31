// bbla_demo.d

import std.stdio;
import bbla;

void main()
{
    writeln("Bare-bones linear algebra demo...");
    writeln("First, some matrix handling.");
    Matrix a = eye(3);
    writeln("a=", a);
    Matrix b = zeros(2,3);
    Matrix c = transpose(b);
    b[1,2] = 99.0;
    c[0,0] = 1.0; c[1,1] = 1.0;
    writeln("b=", b);
    writeln("c=", c);

    Matrix e = new Matrix([1.0, 2.0, 3.0]);
    writeln("e= ", e);
    Matrix e2 = new Matrix([1, 2, 3], "row");
    writeln("e2= ", e2);
    Matrix f = new Matrix([[1.0,2.0,3.0],[4.0,5.0,6.0]]);
    writeln("f= ", f);
    writeln("2*f= ", 2*f);
    writeln("f+2= ", f+2);
    writeln("hstack([f,2*f,3+f])= ", vstack([f,2*f,3+f]));

    Matrix g = dot(f,c);
    writeln("g= ", g);

    writeln("Now, try some linear algebra.");
    writeln("Gauss-Jordan elimination");
    Matrix A = new Matrix([[0.0,  2.0,  0.0,  1.0],
			   [2.0,  2.0,  3.0,  2.0],
			   [4.0, -3.0,  0.0,  1.0],
			   [6.0,  1.0, -6.0, -5.0]]);
    b = new Matrix([[ 0.0],
		    [-2.0],
		    [-7.0],
		    [ 6.0]]);
    writeln("A= ", A, " b= ", b);
    Matrix Ab = hstack([A,b]);
    writeln("before elimination Ab= ", Ab);
    Matrix Aonly = Ab.sliceDup(0, 4, 0, 4);
    Matrix bonly = Ab.sliceDup(0, 4, 4, 5);
    gaussJordanElimination(Ab);
    writeln("after elimination Ab= ", Ab);
    double[] x = Ab.getColumn(4);
    writeln("x= ", x);
    Matrix rhs = dot(Aonly, new Matrix(x));
    writeln("rhs= ", rhs);
    Matrix residual = rhs - bonly;
    writeln("residual= ", residual);

    writeln("LU decomposition");
    c = new Matrix(A);
    auto perm = decomp(c);
    writeln("after decomp, c= ", c);
    writeln("perm= ", perm);
    solve(c, b, perm);
    writeln("solution vector= ", b);

    auto y = inverse(A);
    writeln("dot(A,inverse(A))= ", dot(A,y));

    writeln("Least-squares solve of an over-constrained system.");
    A = new Matrix([[0.0,  2.0,  0.0,  1.0],
		    [2.0,  2.0,  3.0,  2.0],
		    [4.0,  4.0,  6.0,  4.0],
		    [4.0, -3.0,  0.0,  1.0],
		    [4.0, -3.0,  0.0,  1.0],
		    [6.0,  1.0, -6.0, -5.0]]);
    b = new Matrix([[ 0.0],
		    [-2.0],
		    [-4.0],
		    [-7.0],
		    [-7.0],
		    [ 6.0]]);
    writeln("A= ", A, " b= ", b);
    auto xx = lsqsolve(A, b);
    writeln("xx=", xx);

    writeln("Done.");
}

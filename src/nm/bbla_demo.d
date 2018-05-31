// bbla_demo.d

import std.stdio;
import nm.bbla;

void main()
{
    writeln("Bare-bones linear algebra demo...");
    writeln("First, some matrix handling.");
    auto a = eye!double(3);
    writeln("a=", a);
    auto b = zeros!double(2,3);
    auto c = transpose!double(b);
    b[1,2] = 99.0;
    c[0,0] = 1.0; c[1,1] = 1.0;
    writeln("b=", b);
    writeln("c=", c);

    auto e = new Matrix!double([1.0, 2.0, 3.0]);
    writeln("e= ", e);
    auto e2 = new Matrix!double([1, 2, 3], "row");
    writeln("e2= ", e2);
    auto f = new Matrix!double([[1.0,2.0,3.0],[4.0,5.0,6.0]]);
    writeln("f= ", f);
    writeln("2*f= ", 2*f);
    writeln("f+2= ", f+2);
    writeln("vstack([f,2*f,3+f])= ", vstack!double([f,2*f,3+f]));

    auto g = dot!double(f,c);
    writeln("g= ", g);

    writeln("Now, try some linear algebra.");
    writeln("Gauss-Jordan elimination");
    auto A = new Matrix!double([[0.0,  2.0,  0.0,  1.0],
                                [2.0,  2.0,  3.0,  2.0],
                                [4.0, -3.0,  0.0,  1.0],
                                [6.0,  1.0, -6.0, -5.0]]);
    b = new Matrix!double([[ 0.0],
                           [-2.0],
                           [-7.0],
                           [ 6.0]]);
    writeln("A= ", A, " b= ", b);
    auto Ab = hstack!double([A,b]);
    writeln("before elimination Ab= ", Ab);
    auto Aonly = Ab.sliceDup(0, 4, 0, 4);
    auto bonly = Ab.sliceDup(0, 4, 4, 5);
    gaussJordanElimination!double(Ab);
    writeln("after elimination Ab= ", Ab);
    double[] x = Ab.getColumn(4);
    writeln("x= ", x);
    auto rhs = dot!double(Aonly, new Matrix!double(x));
    writeln("rhs= ", rhs);
    auto residual = rhs - bonly;
    writeln("residual= ", residual);

    writeln("LU decomposition");
    c = new Matrix!double(A);
    auto perm = decomp!double(c);
    writeln("after decomp, c= ", c);
    writeln("perm= ", perm);
    solve!double(c, b, perm);
    writeln("solution vector= ", b);

    auto y = inverse!double(A);
    writeln("dot(A,inverse(A))= ", dot!double(A,y));

    writeln("Least-squares solve of an over-constrained system.");
    A = new Matrix!double([[0.0,  2.0,  0.0,  1.0],
                           [2.0,  2.0,  3.0,  2.0],
                           [4.0,  4.0,  6.0,  4.0],
                           [4.0, -3.0,  0.0,  1.0],
                           [4.0, -3.0,  0.0,  1.0],
                           [6.0,  1.0, -6.0, -5.0]]);
    b = new Matrix!double([[ 0.0],
                           [-2.0],
                           [-4.0],
                           [-7.0],
                           [-7.0],
                           [ 6.0]]);
    writeln("A= ", A, " b= ", b);
    auto xx = lsqsolve!double(A, b);
    writeln("xx=", xx);

    writeln("Done.");
}

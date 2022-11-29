// rsla_test.cu
// Test program for Really-simple linear algebra module.
// PJ 2022-10-17, 2022-11-29
// Build with
// $ nvcc --expt-relaxed-constexpr -o rsla_test rsla_test.cu
// Run with
// $ ./rsla_test
//
#include <iostream>
#include "../number.cu"
#include "../rsla.cu"

using namespace std;

int main()
{
    cout << "Test for rsla." << endl;
    cout << "2x2" << endl;
    number A[2][2]{0.0, 2.0, 2.0, 2.0};
    number x[2]{-0.5, 1.0};
    number b[2]; MVMult(A,x, b);
    cout << "A=" << toString(A) << endl;
    cout << "x=" << toString2(x) << endl;
    cout << "b=" << toString2(b) << endl;
    number Ainv[2][2]; int status = MInverse(A, Ainv, 1.0e-9);
    cout << "status=" << status << endl;
    cout << "Ainv=" << toString(Ainv) << endl;
    number y[2]; MVMult(Ainv,b, y);
    cout << "y=" << toString2(y) << endl;
    //
    cout << "and again, with a 3x3 system." << endl;
    number x3[3]{-0.5, 1.0, 1.0/3};
    number A33[3][3]{0.0, 2.0, 0.0,  2.0, 2.0, 3.0,  4.0, -3.0, 0.0};
    number b3[3]; MVMult(A33, x3, b3);
    cout << "A33=" << toString(A33) << endl;
    cout << "x3=" << toString3(x3) << endl;
    cout << "b3=" << toString3(b3) << endl;
    number A33inv[3][3]; status = MInverse(A33, A33inv);
    cout << "A33inv=" << toString(A33) << endl;
    number y3[3]; MVMult(A33inv, b3, y3);
    cout << "y3=" << toString3(y3) << endl;
    //
    cout << "Try Matrix bits, size 3." << endl;
    Matrix<3,1> xx{-0.5, 1.0, 1.0/3};
    Matrix<3,3> AA{0.0, 2.0, 0.0,  2.0, 2.0, 3.0,  4.0, -3.0, 0.0};
    cout << "AA=" << AA << endl;
    cout << "xx=" << xx << endl;
    Matrix<3,1> bb = dot<3,3,1>(AA, xx);
    cout << "bb=" << bb << endl;
    Matrix<3,4> cc = hstack(AA, bb);
    cout << "Augmented Matrix cc=hstack(AA,bb)=" << cc << endl;
    if (gaussJordanElimination(cc)) {
        cout << "Matrix was singular." << endl;
    }
    cout << "After Gauss-Jordan elimination cc=" << cc << endl;
    number yy[3];
    if(cc.getColumn(3, yy)) { cout << "Invalid column." << endl; }
    cout << "yy=["; for (int i=0; i < 3; ++i) cout << yy[i] << ((i < 2) ? "," : "]"); cout << endl;
    //
    cout << "Try Matrix bits, size 4." << endl;
    Matrix<4,1> bbb{0.0, -2.0, -7.0, 6.0};
    Matrix<4,4> AAA {
            0.0,  2.0,  0.0,  1.0,
            2.0,  2.0,  3.0,  2.0,
            4.0, -3.0,  0.0,  1.0,
            6.0,  1.0, -6.0, -5.0
            };
    cout << "AAA=" << AAA << endl;
    cout << "bbb=" << bbb << endl;
    Matrix<4,5> Ab = hstack(AAA, bbb);
    cout << "Augmented Matrix Ab=hstack(AAA,bbb)=" << Ab << endl;
    if (gaussJordanElimination(Ab)) {
        cout << "Matrix was singular." << endl;
    }
    cout << "After Gauss-Jordan elimination Ab=" << Ab << endl;
    Matrix<4,1> xxx;
    if(Ab.getColumn(4, xxx)) { cout << "Invalid column." << endl; }
    cout << "xxx=" << xxx << endl;
    Matrix<4,1> bbb2 = dot<4,4,1>(AAA, xxx);
    cout << "bbb2=" << bbb2 << endl;
    cout << "bbb.approxEqual(bbb2)=" << bbb.approxEqual(bbb2) << endl;
    //
    return 0;
}

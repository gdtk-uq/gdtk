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
using namespace rsla;

int main()
{
    cout << "Test for rsla." << endl;
    cout << "2x2" << endl;
    number A[2][2]{0.0, 2.0, 2.0, 2.0};
    number x[2]{-0.5, 1.0};
    number b[2]; MVMult(A,x, b);
    cout << "A=" << toString<2,2>(A) << endl;
    cout << "x=" << toString<2>(x) << endl;
    cout << "b=" << toString<2>(b) << endl;
    number Ainv[2][2]; int status = MInverse(A, Ainv, 1.0e-9);
    cout << "status=" << status << endl;
    cout << "Ainv=" << toString<2,2>(Ainv) << endl;
    number y[2]; MVMult(Ainv,b, y);
    cout << "y=" << toString<2>(y) << endl;
    //
    cout << "and again, with a 3x3 system." << endl;
    number x3[3]{-0.5, 1.0, 1.0/3};
    number A33[3][3]{0.0, 2.0, 0.0,  2.0, 2.0, 3.0,  4.0, -3.0, 0.0};
    number b3[3]; MVMult(A33, x3, b3);
    cout << "A33=" << toString<3,3>(A33) << endl;
    cout << "x3=" << toString<3>(x3) << endl;
    cout << "b3=" << toString<3>(b3) << endl;
    number A33inv[3][3]; status = MInverse(A33, A33inv);
    cout << "A33inv=" << toString<3,3>(A33) << endl;
    number y3[3]; MVMult(A33inv, b3, y3);
    cout << "y3=" << toString<3>(y3) << endl;
    //
    cout << "Try Matrix bits, size 3." << endl;
    number xx[3][1] {-0.5, 1.0, 1.0/3};
    number AA[3][3] {0.0, 2.0, 0.0,  2.0, 2.0, 3.0,  4.0, -3.0, 0.0};
    cout << "AA=" << toString<3,3>(AA) << endl;
    cout << "xx=" << toString<3,1>(xx) << endl;
    number bb[3][1]; dot<3,3,1>(AA, xx, bb);
    cout << "bb=" << toString<3,1>(bb) << endl;
    number cc[3][4]; hstack<3,3,1>(AA, bb, cc);
    cout << "Augmented Matrix cc=hstack(AA,bb)=" << toString<3,4>(cc) << endl;
    if (gaussJordanElimination<3,4>(cc)) {
        cout << "Matrix was singular." << endl;
    }
    cout << "After Gauss-Jordan elimination cc=" << toString<3,4>(cc) << endl;
    number yy[3];
    if(getColumn<3,4>(cc, 3, yy)) { cout << "Invalid column." << endl; }
    cout << "yy=" << toString<3>(yy) << endl;
    //
    cout << "Try Matrix bits, size 4." << endl;
    number bbb[4][1] {0.0, -2.0, -7.0, 6.0};
    number AAA[4][4] {
            0.0,  2.0,  0.0,  1.0,
            2.0,  2.0,  3.0,  2.0,
            4.0, -3.0,  0.0,  1.0,
            6.0,  1.0, -6.0, -5.0
            };
    cout << "AAA=" << toString<4,4>(AAA) << endl;
    cout << "bbb=" << toString<4,1>(bbb) << endl;
    number Ab[4][5]; hstack<4,4,1>(AAA, bbb, Ab);
    cout << "Augmented Matrix Ab=hstack(AAA,bbb)=" << toString<4,5>(Ab) << endl;
    if (gaussJordanElimination<4,5>(Ab)) {
        cout << "Matrix was singular." << endl;
    }
    cout << "After Gauss-Jordan elimination Ab=" << toString<4,5>(Ab) << endl;
    number xxx[4][1];
    if(getColumn<4,5>(Ab, 4, xxx)) { cout << "Invalid column." << endl; }
    cout << "xxx=" << toString<4,1>(xxx) << endl;
    number bbb2[4][1]; dot<4,4,1>(AAA, xxx, bbb2);
    cout << "bbb2=" << toString<4,1>(bbb2) << endl;
    cout << "approxEqual(bbb, bbb2)=" << approxEqual<4,1>(bbb, bbb2) << endl;
    //
    return 0;
}

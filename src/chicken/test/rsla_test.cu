// rsla_test.cu
// Test program for Really-simple linear algebra module.
// PJ 2022-10-17
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
    return 0;
}

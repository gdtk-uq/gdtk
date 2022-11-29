// vector3_test.cu
// Test program for Vector3 module.
// PJ 2022-11-19
// Build with
// $ nvcc --expt-relaxed-constexpr -o vector3_test vector3_test.cu
// Run with
// $ ./vector3_test
//
#include <iostream>
#include "../number.cu"
#include "../vector3.cu"

using namespace std;

int main()
{
    cout << "Test for Vector3." << endl;
    Vector3 a = Vector3{1.0, 2.0, 3.0};
    Vector3 x = Vector3{1.0, 0.0, 0.0};
    Vector3 y = Vector3{0.0, 1.0, 0.0};
    Vector3 z = Vector3{0.0, 0.0, 1.0};
    //
    Vector3 n = y; Vector3 t1 = z; Vector3 t2 = cross(n, t1);
    cout << "t2.approxEqual(x)=" << t2.approxEqual(x) << endl;
    cout << "Original a=" << a << endl;
    a.normalize();
    cout << "Normalized a=" << a << endl;
    a.transform_to_local_frame(n, t1, t2);
    cout << "Local frame" << endl
         << "    n=" << n << endl
         << "   t1=" << t1 << endl
         << "   t2=" << t2 << endl;
    cout << "In local frame  a=" << a << endl;
    a.transform_to_global_frame(n, t1, t2);
    cout << "In global frame a=" << a << endl;
    //
    return 0;
}

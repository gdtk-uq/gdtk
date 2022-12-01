// spline_test.cu
// Test program for the CubicSpline.
//
// PJ 2022-11-30

#include <iostream>
#include "../number.cu"
#include "../spline.cu"

using namespace std;

number runge(number x) { return 1.0/(1.0 + 25* x * x); }

int main()
{
    cout << "Test for CubicSpline." << endl;
    cout << "Simple initialization example." << endl;
    number x[5]{0.0, 1.0, 2.0, 3.0, 4.0};
    number y[5]{0.0, 1.0, 0.0, 1.0, 0.0};
    CubicSpline<4> s0; s0.set(x, y);
    cout << "s0=" << s0.toString() << endl;

    cout << "Runge function example." << endl;
    constexpr int ns = 21;
    number x0 = -1.0;
    number x1 = 1.0;
    number dx = (x1-x0)/(ns-1);
    number x_sample[ns], y_sample[ns];
    for (int i=0; i < ns; ++i) {
        number xx = x0 + dx*i;
        x_sample[i] = xx;
        y_sample[i] = runge(xx);
    }
    CubicSpline<ns-1> s; s.set(x_sample, y_sample);
    cout << "s=" << s.toString() << endl;
    constexpr int ms = 100;
    dx = (x1-x0)/(ms-1);
    number max_dy = 0.0;
    for (int i=0; i < ms; ++i) {
        number xx = x0 + dx*i;
        number y_runge = runge(xx);
        number y_spline = s.eval(xx);
        number dy = y_spline - y_runge;
        max_dy = (abs(dy) > max_dy) ? abs(dy) : max_dy;
        // cout << xx << " " << y_runge << " " << y_spline << " " << dy << endl;
        // assert(fabs(dy) < 0.02, failedUnitTest());
    }
    cout << "max_dy=" << max_dy << endl;
    return 0;
}

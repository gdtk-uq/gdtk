// ridder_demo.d
// Solve a nonlinear equation f(x)=0 using the method of Ridder.
// Peter J. 
// demo code 13-Jun-2014, adapted from the C++ version

import std.stdio;
import std.math;
import ridder;

double test_fun_1(double x) {
    return pow(x,3) + pow(x,2) - 3*x - 3;
}

double test_fun_2(double x, double a) {
    return a*x + sin(x) - exp(x);
}

void main() {
    writeln("Begin demonstration of using Ridder's function solver...");
    writeln();
    writeln("Example 1 is from Gerald and Wheatley, p. 45");
    writeln("Solve f(x) = x^3 + x^2 - 3x -3 = 0 with initial ",
            "guesses of x0 = 1 and x1 = 2.");
    writeln("Final result x = ", solve!(test_fun_1)(1, 2));
    writeln("Gerald and Wheatley report x = 1.732051");
    writeln();
    //
    writeln("Example 2 is from Gerald and Wheatley, p.45 also, ", 
            "but using a delegate for closure.");
    writeln("Solve f(x) = 3*x + sin(x) - e^x = 0 with initial ",
            "guesses of x0 = 0 and x1 = 1.");
    double my_a = 3.0;
    auto test_fun_3 = delegate (double x) { return test_fun_2(x, my_a); }; 
    writeln("Final result x = ", solve!test_fun_3(0, 1));
    writeln("Gerald and Wheatley report x = 0.3604217");
    writeln();
    //
    writeln("Bracket a root of the second example.");
    double x1 = 0.4;
    double x2 = 0.5;
    int result_flag = bracket!test_fun_3(x1, x2);
    writeln("result_flag=", result_flag, " x1=", x1, " x2=", x2);
    writeln("Done.");
}

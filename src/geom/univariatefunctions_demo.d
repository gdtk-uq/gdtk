// univariatefunctions_demo.d

import std.stdio;
import univariatefunctions;

void main()
{
    writeln("Begin demonstration of the clustering functions.");
    auto cf = new RobertsFunction(false, true, 1.1);
    auto cf2 = new LinearFunction(1.0, 0.0);
    foreach(i; 0 .. 11) {
	double x = 0.1*i;
	writeln("x= ", x, " roberts_y= ", cf(x), " linear_y= ", cf2(x));
    }
} // end main

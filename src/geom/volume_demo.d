/**
 * volume_demo.d Demonstrate some of the behaviour of constructing volumes.
 *
 * Author: Peter J.
 * Version: 2015-04-07
 */

import std.stdio;
import geom;
import gpath;
import surface;
import volume;

void main()
{
    writeln("Begin demonstration of the ParametricVolume elements.");
    Vector3[8] p;
    p[0] = Vector3(0.0, 0.1, 0.0);
    p[1] = Vector3(1.0, 0.1, 0.0);
    p[2] = Vector3(1.0, 1.1, 0.0);
    p[3] = Vector3(0.0, 1.1, 0.0);

    p[4] = Vector3(0.0, 0.1, 3.0);
    p[5] = Vector3(1.0, 0.1, 3.0);
    p[6] = Vector3(1.0, 1.1, 3.0);
    p[7] = Vector3(0.0, 1.1, 3.0);

    writeln("TFIVolume demo");
    auto my_vol = new TFIVolume(p);
    writeln("my_vol= ", my_vol);
    auto c = my_vol(0.1, 0.1, 0.1);
    writeln("my_vol(0.1, 0.1, 0.1)= ", c);

    writeln("SubRangedVolume demo");
    auto srv = new SubRangedVolume(my_vol, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5);
    writeln("srv(0.2,0.2,0.2)= ", srv(0.2,0.2,0.2));
    writeln("Done volume_demo.");
}

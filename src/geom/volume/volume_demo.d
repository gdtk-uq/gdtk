/**
 * volume_demo.d Demonstrate some of the behaviour of constructing volumes.
 *
 * Author: Peter J.
 * Version: 2015-04-07
 */

import std.stdio;
import geom;

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

    writeln("SweptSurfaceVolume demo");
    auto face0123 = new CoonsPatch(p[0], p[1], p[2], p[3]);
    auto edge04 = new Line(p[0], p[4]);
    auto ssv = new SweptSurfaceVolume(face0123, edge04);
    writeln("ssv=", ssv);
    writeln("ssv(0.1, 0.1, 0.1)= ", ssv(0.1, 0.1, 0.1));

    writeln("TwoSurfaceVolume demo");
    face0123 = new CoonsPatch(p[0], p[1], p[2], p[3]);
    auto face4567 = new CoonsPatch(p[4], p[5], p[6], p[7]);
    auto tsv = new TwoSurfaceVolume(face0123, face4567);
    writeln("tsv=", tsv);
    writeln("tsv(0.1, 0.1, 0.1)= ", tsv(0.1, 0.1, 0.1));

    writeln("Done volume_demo.");
}

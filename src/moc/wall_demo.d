/*
 * file: moc_wall_demo.d
 * location: /home/momar/Documents/Thesis/IMOC_D
 * author: Momar Hughes
 * date: 30 Apr 2015
 * description: Demo of moc_wall.d
 */

// dmd wall_demo.d gasdynamic.d kernel.d wall.d geom.d gpath.d linesearch.d bbla.d

import std.stdio, std.math;
import wall;
import geom;
import gpath;

void main(){
    writeln("Begin moc_wall_demo.d");

    writeln("Test: WallFromPoints, DeleteWall, WallIsPresent...");
    auto points = [Vector3(0.0,0.0),Vector3(10.0,5.0)];
    WallFromPoints(0,"Line",points);
    writeln("... wall made");
    writefln("Wall 0 deleted: %d, (1 for success)",WallIsPresent(0));
    writeln();
    /*writeln("Test: WallFromPath...");
      Line path = new Line(Vector3(0.0,0.0),Vector3(10.0,5.0));
      WallFromPath(0,path);
      writefln("Wall 0 is present: %d, (1 for success)",WallIsPresent(0));*/
    writefln("Wall 0 has %d points: (%s),(%s)",WallGetNumberOfPoints(0),WallGetPoint(0,0),WallGetPoint(0,1));
    writeln("...done");
    writeln();
    writeln("Test: WallPos, WallSlope");
    double t = 0.75;
    writefln("Coordinate of wall for t=%g: %s",t,WallPos(0,t));
    writefln("the slope here is %g",WallSlope(0,t));
    writeln("...done");
    writeln();
        
    writeln("Test: WallFindT...");
    auto a = Vector3(0.0,3.0);
    double cosine = 0.707, sine = -0.707;
    t = WallFindT(0,a,cosine,sine);
    writefln("Point: %s, cos = %g, sin = %g",a,cosine,sine);
    writefln("The intersection is at %s, t=%g",WallPos(0,t),t);
    writeln("...done");
    writeln();
        
    writeln("Test: SaveWall...");
    SaveWall(0,"moc_wall_demo.txt");
    writeln("... wall 0 saved.");
    DeleteWall(0);
        
    writeln("Test: LoadWall...");
    LoadWall(0,"moc_wall_demo.txt");
    writeln("... wall 0 loaded.");
    writefln("Wall 0 is present: %d, (1 for success)",WallIsPresent(0));
    writefln("Wall 0 has %d points: (%s),(%s)",WallGetNumberOfPoints(0),WallGetPoint(0,0),WallGetPoint(0,1));
    writeln();
    writeln("End moc_wall_demo.d");
} // end main()
        

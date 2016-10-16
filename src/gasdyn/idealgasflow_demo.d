/*
 * file: idealgasflow_demo.d
 * author: Momar Hughes
 * date: 5 Apr 2015 initial code
 */

import std.stdio, std.math;
import idealgasflow;

void main(){
    writeln("Begin gasdynamic demo...");
    double M = 2.4, g = 1.4;
    writefln("Use M=%g and gamma=%g",M,g);
	
    writeln("Test: T0_T...");
    writeln("...returns T0/T = ",T0_T(M,g));
    writeln("Literature value = 2.152");
	
    writeln("Test: p0_p...");
    writeln("...returns p0/p = ",p0_p(M,g));
    writeln("Literature value = 14.620");
	
    writeln("Test: NuFromM...");
    writeln("...returns nu = ",PM1(M,g));
    writeln("Literature value = 0.6413");
    // [TODO] These tests for subsonic values should go into the unittest.
    writeln("...provided M=0.8:");
    try {
	PM1(0.8,g);
    } catch (Error e) {
	writeln("Caught subsonic value.");
    }
	
    double nu = 0.6413479572;
    writeln("Test: MFromNu...");
    writeln("...returns M = ",PM2(nu,g));
    writeln("Literature value = 2.4");
    writeln("...provided nu=-0.5:");
    try {
	PM2(-0.5,g);
    } catch (Error e) {
	writeln("Caught subsonic value.");
    }
	
    writeln("Test: MachAngle...");
    writeln("...returns mu = ",MachAngle(M));
    writeln("Literature value = 0.430");
    writeln("...provided M=0.8:");
    try {
	MachAngle(0.8);
    } catch (Error e) {
	writeln("Caught subsonic value.");
    }
	
    writeln("Finished gasdynamic demo.");
} // end main()

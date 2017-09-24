// globalconfig_demo.d

import std.stdio;
import geom;
import globalconfig;

void main()
{
    writeln("test globalconfig module");
    GlobalConfig.nFluidBlocks = 21;
    GlobalConfig.turbulent_zones ~= new BlockZone(Vector3(0.0,0.0,0.0), Vector3(1.0,1.0,1.0));
    GlobalConfig.turbulent_zones ~= new BlockZone(Vector3(0.0,0.0,0.0), Vector3(1.0,1.0,1.0));
    writeln("dimensions= ", GlobalConfig.dimensions);
    writeln("nFluidBlocks= ", GlobalConfig.nFluidBlocks);
    writeln("turbulent_zones= ", GlobalConfig.turbulent_zones);
    writeln("done.");
}


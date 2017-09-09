/**
 * solidfvinterface.d
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-22-04: initial cut
 */

module solidfvinterface;

import geom;
import solidprops;
import ssolidblock;
import globaldata;
import globalconfig;
import solidfvcell;

class SolidFVInterface {
public:
    size_t id;
    // Geometry
    Vector3 pos;
    double Ybar;
    double length;
    double area;
    Vector3 n;
    Vector3 t1;
    Vector3 t2;
    // Material propertis
    SolidProps sp;
    // State
    double T;
    double e;
    double flux;
}

void initPropertiesAtSolidInterfaces(SSolidBlock[] sblks)
{
    if (GlobalConfig.dimensions == 3 ) {
	throw new Error("initPropertiesAtSolidInterfaces not implemented for 3D calculations.");
    }

    SolidFVCell Lft, Rght;
    SolidFVInterface IFace;
    size_t k = 0;

    foreach (blk; sblks) {
	// Do i-facing interfaces
	for (size_t j = blk.jmin; j <= blk.jmax; ++j) {
	    for (size_t i = blk.imin; i <= blk.imax+1; ++i) {
		Lft = blk.getCell(i-1,j,k);
		Rght = blk.getCell(i,j,k);
		if ( i == blk.imin ) {
		    // There is no Lft cell, so
		    Lft = Rght;
		}
		if ( i == blk.imax+1 ) {
		    // There is no right cell, so
		    Rght = Lft;
		}
		IFace = blk.getIfi(i,j,k);
		IFace.sp = new SolidProps(Lft.sp, Rght.sp, 0.5, 0.5);
	    }
	}
	// Do j-facing interfaces
	for (size_t j = blk.jmin; j <= blk.jmax+1; ++j) {
	    for (size_t i = blk.imin; i <= blk.imax; ++i) {
		Lft = blk.getCell(i,j-1,k);
		Rght = blk.getCell(i,j,k);
		if ( j == blk.jmin ) {
		    // There is no Lft cell, so
		    Lft = Rght;
		}
		if ( j == blk.jmax+1 ) {
		    // There is no right cell, so
		    Rght = Lft;
		}
		IFace = blk.getIfj(i,j,k);
		IFace.sp = new SolidProps(Lft.sp, Rght.sp, 0.5, 0.5);
	    }
	}
    }
}



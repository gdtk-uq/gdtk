/**
 * solidfvinterface.d
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-22-04: initial cut
 */

module solidfvinterface;

import nm.complex;
import nm.number;

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
    number Ybar;
    number length;
    number area;
    Vector3 n;
    Vector3 t1;
    Vector3 t2;
    // Material propertis
    SolidProps sp;
    // State
    number T;
    number e;
    number flux;
}

void initPropertiesAtSolidInterfaces(SSolidBlock[] sblks)
{
    //if (GlobalConfig.dimensions == 3 ) {
    //    throw new Error("initPropertiesAtSolidInterfaces not implemented for 3D calculations.");
    //}

    SolidFVCell Lft, Rght;
    SolidFVInterface IFace;

    /*
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
                IFace.sp = SolidProps(Lft.sp, Rght.sp, 0.5, 0.5);
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
                IFace.sp = SolidProps(Lft.sp, Rght.sp, 0.5, 0.5);
            }
        }
    }
    */
    foreach (blk; sblks) {
        // ifi interfaces are west interfaces, with their unit normal pointing east.
        for (size_t k = blk.kmin; k <= blk.kmax; ++k) {
            for (size_t j = blk.jmin; j <= blk.jmax; ++j) {
                for (size_t i = blk.imin; i <= blk.imax+1; ++i) {
                    IFace = blk.getIfi(i,j,k);
                    Lft = (blk.myConfig.dimensions == 3) ? blk.getCell(i-1,j,k) : blk.getCell(i-1,j);
                    Rght = (blk.myConfig.dimensions == 3) ? blk.getCell(i,j,k) : blk.getCell(i,j);
                    if (i == blk.imin) {
                        // There is no Lft cell, so
                        Lft = Rght;
                    }
                    if (i == blk.imax+1) {
                        // There is no Rght cell, so
                        Rght = Lft;
                    }
                    IFace.sp = SolidProps(Lft.sp, Rght.sp, 0.5, 0.5);
                } // i loop
            } // j loop
        } // for k
        // ifj interfaces are south interfaces, with their unit normal pointing north.
        for (size_t k = blk.kmin; k <= blk.kmax; ++k) {
            for (size_t i = blk.imin; i <= blk.imax; ++i) {
                for (size_t j = blk.jmin; j <= blk.jmax+1; ++j) {
                    IFace = blk.getIfj(i,j,k);
                    Lft = (blk.myConfig.dimensions == 3) ? blk.getCell(i,j-1,k) : blk.getCell(i,j-1);
                    Rght = (blk.myConfig.dimensions == 3) ? blk.getCell(i,j,k) : blk.getCell(i,j);
                    if (j == blk.jmin) {
                        // There is no Lft cell, so
                        Lft = Rght;
                    }
                    if (j == blk.jmax+1) {
                        // There is no Rght cell, so
                        Rght = Lft;
                    }
                    IFace.sp = SolidProps(Lft.sp, Rght.sp, 0.5, 0.5);
                } // j loop
            } // i loop
        } // for k
        if (blk.myConfig.dimensions == 3) {
            // ifk interfaces are bottom interfaces, with unit normal pointing to top.
            for (size_t i = blk.imin; i <= blk.imax; ++i) {
                for (size_t j = blk.jmin; j <= blk.jmax; ++j) {
                    for (size_t k = blk.kmin; k <= blk.kmax+1; ++k) {
                        IFace = blk.getIfk(i,j,k);
                        Lft = blk.getCell(i,j,k-1);
                        Rght = blk.getCell(i,j,k);
                        if (k == blk.kmin) {
                            // There is no Lft cell, so
                            Lft = Rght;
                        }
                        if (k == blk.kmax+1) {
                            // There is no Rght cell, so
                            Rght = Lft;
                        }
                        IFace.sp = SolidProps(Lft.sp, Rght.sp, 0.5, 0.5);
                    } // for k
                } // j loop
            } // i loop
        }
    }
}

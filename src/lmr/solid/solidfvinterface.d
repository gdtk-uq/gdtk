/**
 * solidfvinterface.d
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-22-04: initial cut
 */

module solidfvinterface;

import std.math;

import ntypes.complex;
import nm.number;

import geom;
import ssolidblock;
import globaldata;
import globalconfig;
import solidfvcell;
import lmr.solid.solidstate;

class SolidFVInterface {
public:
    size_t id;
    bool is_on_boundary = false;  // by default, assume not on boundary
    size_t bc_id;   // if the face is on a block boundary, which one
    size_t bc_idx;  // if the face is on a block boundary, store index into the array of faces attached to bc
    // Geometry
    Vector3 pos;
    number Ybar;
    number length;
    number area;
    Vector3 n;
    Vector3 t1;
    Vector3 t2;
    // Inteface state
    SolidState ss;
    number T;
    number e;
    number flux;
    // Temperature gradient
    number dTdx;
    number dTdy;
    number dTdz;
    // Cells attached to interface
    SolidFVCell cellLeft;
    SolidFVCell cellRight;

    @nogc
    void averageTemperature() {
        // set the inteface temperature to the arithmetic average of the two adjoining cells.
        if (!cellLeft.is_ghost && !cellRight.is_ghost) { T = 0.5*(cellLeft.T + cellRight.T); }
    }

    @nogc
    void averageProperties() {
        if (!cellLeft.is_ghost && !cellRight.is_ghost) {
            harmonicAverage(cellLeft.ss, cellRight.ss, ss);
        }
    }

    @nogc
    void averageTGradient() {
        // set the interface gradients to the arithmetic average of the two adjoining cells
        SolidFVCell cL0 = cellLeft;
        SolidFVCell cR0 = cellRight;
        if (cL0.is_ghost) {
            dTdx = cR0.dTdx;
            dTdy = cR0.dTdy;
            dTdz = cR0.dTdz;
        } else if (cR0.is_ghost) {
            dTdx = cL0.dTdx;
            dTdy = cL0.dTdy;
            dTdz = cL0.dTdz;
        } else {
            dTdx = 0.5*(cL0.dTdx+cR0.dTdx);
            dTdy = 0.5*(cL0.dTdy+cR0.dTdy);
            dTdz = 0.5*(cL0.dTdz+cR0.dTdz);
        }
    }

    @nogc
    void augmentTGradient() {
        // augment the arithmetic average to help prevent odd-even decoupling.
        SolidFVCell cL0 = cellLeft;
        SolidFVCell cR0 = cellRight;
        // interface normal
        number nx = n.x;
        number ny = n.y;
        number nz = n.z;
        // vector from left-cell-centre to right-cell-centre
        number ex = cR0.pos.x - cL0.pos.x;
        number ey = cR0.pos.y - cL0.pos.y;
        number ez = cR0.pos.z - cL0.pos.z;
        // ehat
        number emag = sqrt(ex*ex + ey*ey + ez*ez);
        number ehatx = ex/emag;
        number ehaty = ey/emag;
        number ehatz = ez/emag;
        // ndotehat
        number ndotehat = nx*ehatx + ny*ehaty + nz*ehatz;
        number avgdotehat = 0.5*(cL0.dTdx+cR0.dTdx)*ehatx +
            0.5*(cL0.dTdy+cR0.dTdy)*ehaty +
            0.5*(cL0.dTdz+cR0.dTdz)*ehatz;
        number jump = avgdotehat - (cR0.T - cL0.T)/emag;
        if (!cL0.is_ghost && !cR0.is_ghost) {
            dTdx -= jump*(nx/ndotehat);
            dTdy -= jump*(ny/ndotehat);
            dTdz -= jump*(nz/ndotehat);
        }
    }

    @nogc
    void computeFlux(size_t dim, bool solid_has_isotropic_properties) {
        // Evaluate the flux of energy passing through each interface.
        // Note: this function assumes that the temperature gradients
        //       have already been approximated at the interfaces.
        number qx, qy, qz;
        qx = -ss.k * dTdx;
        qy = -ss.k * dTdy;
        qz = -ss.k * dTdz;
        flux = qx * n.x + qy * n.y + qz * n.z;
    }

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
                    harmonicAverage(Lft.ss, Rght.ss, IFace.ss);
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
                    harmonicAverage(Lft.ss, Rght.ss, IFace.ss);
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
                        harmonicAverage(Lft.ss, Rght.ss, IFace.ss);
                    } // for k
                } // j loop
            } // i loop
        }
    }
}

/**
 * gas_solid_interface.d
 *
 * Authors: RG & PJ
 * Date: 2015-05-07
 *
 * History: This source code is a variation of that
 *          found ind eilmer3/source/conj-ht-interface.cxx.
 *          That code was first worked on by Delphine Francois
 *          and Justin Beri.
 */

module gas_solid_interface;

import std.stdio;
import std.math;
import nm.complex;
import nm.number;

import fvcell;
import fvinterface;
import solidfvcell;
import solidfvinterface;

void computeFluxesAndTemperatures(int ftl, FVCell[] gasCells, FVInterface[] gasIFaces,
                                  SolidFVCell[] solidCells, SolidFVInterface[] solidIFaces)
{
    number dxG, dyG, dnG, dxS, dyS, dnS;
    number kG_dnG, kS_dnS, cosA, cosB;
    number T, q;

    foreach ( i; 0 .. gasCells.length ) {
        cosA = gasIFaces[i].n.x;
        cosB = gasIFaces[i].n.y;

        dxG = gasIFaces[i].pos.x - gasCells[i].pos[0].x;
        dyG = gasIFaces[i].pos.y - gasCells[i].pos[0].y;
        dnG = fabs(cosA*dxG + cosB*dyG);

        dxS = solidIFaces[i].pos.x - solidCells[i].pos.x;
        dyS = solidIFaces[i].pos.y - solidCells[i].pos.y;
        dnS = fabs(cosA*dxS + cosB*dyS);

        kG_dnG = gasCells[i].fs.gas.k / dnG;
        kS_dnS = solidCells[i].sp.k / dnS;

        T = (gasCells[i].fs.gas.T*kG_dnG + solidCells[i].T*kS_dnS) / (kG_dnG + kS_dnS);
        q = -kG_dnG * (T - gasCells[i].fs.gas.T);
        // Finally update properties in interfaces
        gasIFaces[i].fs.gas.T = T;
        gasIFaces[i].F.total_energy = q; // CHECK ME: might only work for
                                         // NORTH-SOUTH orientation.
                                         // Need to think about sign.
        // Note that we are not propagating complex numbers into solid domain.
        solidIFaces[i].T = T.re;
        solidIFaces[i].flux = q.re;
    }
}

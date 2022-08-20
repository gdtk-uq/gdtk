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
import nm.bbla;

import fvcell;
import fvinterface;
import solidfvcell;
import solidfvinterface;

@nogc
void computeFluxesAndTemperatures(int ftl, FVCell[] gasCells, FVInterface[] gasIFaces,
                                  SolidFVCell[] solidCells, SolidFVInterface[] solidIFaces)
{
    number dxG, dyG, dzG, dnG, dxS, dyS, dzS, dnS;
    number kG_dnG, kS_dnS, cosA, cosB, cosC;
    number T, q;
    auto cqi = gasCells[0].myConfig.cqi;
    //
    foreach ( i; 0 .. gasCells.length ) {
        cosA = gasIFaces[i].n.x;
        cosB = gasIFaces[i].n.y;
        cosC = gasIFaces[i].n.z;
        //
        dxG = gasIFaces[i].pos.x - gasCells[i].pos[0].x;
        dyG = gasIFaces[i].pos.y - gasCells[i].pos[0].y;
        dzG = gasIFaces[i].pos.z - gasCells[i].pos[0].z;
        dnG = fabs(cosA*dxG + cosB*dyG + cosC*dzG);
        //
        dxS = solidIFaces[i].pos.x - solidCells[i].pos.x;
        dyS = solidIFaces[i].pos.y - solidCells[i].pos.y;
        dzS = solidIFaces[i].pos.z - solidCells[i].pos.z;
        dnS = fabs(cosA*dxS + cosB*dyS + cosC*dzS);
        //
        kG_dnG = gasCells[i].fs.gas.k / dnG;
        kS_dnS = solidCells[i].sp.k / dnS;
        //
        T = (gasCells[i].fs.gas.T*kG_dnG + solidCells[i].T*kS_dnS) / (kG_dnG + kS_dnS);
        q = -kG_dnG * (T - gasCells[i].fs.gas.T);
        // Finally update properties in interfaces
        gasIFaces[i].fs.gas.T = T;
        gasIFaces[i].F[cqi.totEnergy] = q; // CHECK ME: might only work for
                                               // NORTH-SOUTH orientation.
                                               // Need to think about sign.
        solidIFaces[i].T = T;
        solidIFaces[i].flux = q;
    }
}

// not @nogc because of LUDecomp
void computeFluxesAndTemperatures2(int ftl, FVCell[] gasCells, FVInterface[] gasIFaces, 
                                   SolidFVCell[] solidCells, SolidFVInterface[] solidIFaces,
                                   number[] T, number[] B, Matrix!number A, int[] pivot)
{
    // 1. Assemble matrix of coefficients
    A.zeros;
    // 1a. First row in matrix (special)
    auto dyG = gasIFaces[0].pos.y - gasCells[0].pos[0].y;
    auto dyS = solidCells[0].pos.y - solidIFaces[0].pos.y;
    auto dx = 2*(solidCells[1].pos.x - solidCells[0].pos.x);
    auto kG = gasCells[0].fs.gas.k;
    auto k22 = solidCells[0].sp.k22;
    auto k21 = solidCells[0].sp.k21;
    auto TG = gasCells[0].fs.gas.T;
    auto TS = solidCells[0].T;
    A[0,0] = kG/dyG + k22/dyS + k21/dx;
    A[0,1] = -k21/dx;
    B[0] = kG*TG/dyG + k22*TS/dyS;
    
    // 1b. Internal rows
    foreach (i; 1 .. gasCells.length-1) {
        dyG = gasIFaces[i].pos.y - gasCells[i].pos[0].y;
        dyS = solidCells[i].pos.y - solidIFaces[i].pos.y;
        dx = 2*(solidCells[i].pos.x - solidCells[i-1].pos.x);
        kG = gasCells[i].fs.gas.k;
        k22 = solidCells[i].sp.k22;
        k21 = solidCells[i].sp.k21;
        TG = gasCells[i].fs.gas.T;
        TS = solidCells[i].T;
        A[i,i-1] = k21/dx;
        A[i,i] = kG/dyG + k22/dyS;
        A[i,i+1] = -k21/dx;
        B[i] = kG*TG/dyG + k22*TS/dyS;
    }

    // 1c. Last row in matrix (special)
    auto nm1 = gasCells.length-1;
    dyG = gasIFaces[nm1].pos.y - gasCells[nm1].pos[0].y;
    dyS = solidCells[nm1].pos.y - solidIFaces[nm1].pos.y;
    dx = 2*(solidCells[nm1].pos.x - solidCells[nm1-1].pos.x);
    kG = gasCells[nm1].fs.gas.k;
    k22 = solidCells[nm1].sp.k22;
    k21 = solidCells[nm1].sp.k21;
    TG = gasCells[nm1].fs.gas.T;
    TS = solidCells[nm1].T;
    A[nm1,nm1-1] = k21/dx;
    A[nm1,nm1] = kG/dyG + k22/dyS - k21/dx;
    B[nm1] = kG*TG/dyG + k22*TS/dyS;

    // 2. Solve for temperatures.
    LUDecomp!number(A, pivot);
    LUSolve!number(A, pivot, B, T);
    // 3. Place temperatures and fluxes in interfaces
    auto cqi = gasCells[0].myConfig.cqi;
    foreach (i; 0 .. gasCells.length) {
        kG = gasCells[i].fs.gas.k;
        dyG = gasIFaces[i].pos.y - gasCells[i].pos[0].y;
	auto kG_dyG = kG/dyG;
        auto q = -kG_dyG*(T[i] - gasCells[i].fs.gas.T);
        gasIFaces[i].fs.gas.T = T[i];
        gasIFaces[i].F[cqi.totEnergy] = q;
        solidIFaces[i].T = T[i];
        solidIFaces[i].flux = q;
    }
}

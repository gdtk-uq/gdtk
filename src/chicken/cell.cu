// cell.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef CELL_INCLUDED
#define CELL_INCLUDED

#include <string>
#include <sstream>
#include <stdexcept>

#include "number.cu"
#include "vector3.cu"
#include "gas.cu"
#include "flow.cu"
#include "vertex.cu"
#include "face.cu"

using namespace std;

namespace Face {
    // Symbolic names for the faces of the cell and of the block.
    constexpr int iminus = 0;
    constexpr int iplus = 1;
    constexpr int jminus = 2;
    constexpr int jplus = 3;
    constexpr int kminus = 4;
    constexpr int kplus = 5;

    array<string,6> names {"iminus", "iplus", "jminus", "jplus", "kminus", "kplus"};
};

int Face_indx_from_name(string name)
{
    if (name == "iminus") return Face::iminus;
    if (name == "iplus") return Face::iplus;
    if (name == "jminus") return Face::jminus;
    if (name == "jplus") return Face::jplus;
    if (name == "kminus") return Face::kminus;
    if (name == "kplus") return Face::kplus;
    throw runtime_error("Invalid face name: " + name);
}


namespace IOvar {
    // Following the new IO model for Eilmer, we set up the accessor functions
    // for the flow data that is held in the flow data files.
    // These accessor functions are associated with the Cell structure.

    // Keep the following list consistent with the GlobalConfig.iovar_names list
    // in chkn_prep.py and with the symbolic constants just below.
    vector<string> names {"pos.x", "pos.y", "pos.z", "vol",
                              "p", "T", "rho", "e", "a",
                              "vel.x", "vel.y", "vel.z"};

    // We will use these symbols to select the varaible of interest.
    constexpr int posx = 0;
    constexpr int posy = posx + 1;
    constexpr int posz = posy + 1;
    constexpr int vol = posz + 1;
    constexpr int p = vol + 1;
    constexpr int T = p + 1;
    constexpr int rho = T + 1;
    constexpr int e = rho + 1;
    constexpr int a = e + 1;
    constexpr int velx = a + 1;
    constexpr int vely = velx + 1;
    constexpr int velz = vely + 1;
    constexpr int n = velz + 1; // number of symbols that point to the flow variables
}


struct FVCell {
    Vector3 pos; // position of centroid
    number volume;
    number iLength, jLength, kLength; // These lengths are used in the interpolation fns.
    FlowState fs;
    // We will keep connections to the pieces compising the cell as indices
    // into the block's arrays.
    // Although we probably don't need build and keep this data for the structured grid,
    // it simplifies some of the geometry and update code and may ease the use of
    // unstructured grids at a later date.
    array<int,8> vtx{0, 0, 0, 0, 0, 0, 0, 0};
    array<int,6> face{0, 0, 0, 0, 0, 0};

    string toString() {
        ostringstream repr;
        repr << "Cell(pos=" << pos.toString() << ", volume=" << volume;
        repr << ", iLength=" << iLength << ", jLength=" << jLength << ", kLength=" << kLength;
        repr << ", fs=" << fs.toString();
        repr << ", vtx=["; for(auto v : vtx) repr << v << ","; repr << "]";
        repr << ", face=["; for(auto v : face) repr << v << ","; repr << "]";
        repr << ")";
        return repr.str();
    }

    void iovar_set(int i, number val)
    {
        switch (i) {
        case IOvar::posx: pos.x = val; break;
        case IOvar::posy: pos.y = val; break;
        case IOvar::posz: pos.z = val; break;
        case IOvar::vol: volume = val; break;
        case IOvar::p: fs.gas.p = val; break;
        case IOvar::T: fs.gas.T = val; break;
        case IOvar::rho: fs.gas.rho = val; break;
        case IOvar::e: fs.gas.e = val; break;
        case IOvar::a: fs.gas.a = val; break;
        case IOvar::velx: fs.vel.x = val; break;
        case IOvar::vely: fs.vel.y = val; break;
        case IOvar::velz: fs.vel.z = val; break;
        default:
            throw runtime_error("Invalid selection for IOvar: "+to_string(i));
        }
    }

    number iovar_get(int i)
    {
        switch (i) {
        case IOvar::posx: return pos.x;
        case IOvar::posy: return pos.y;
        case IOvar::posz: return pos.z;
        case IOvar::vol: return volume;
        case IOvar::p: return fs.gas.p;
        case IOvar::T: return fs.gas.T;
        case IOvar::rho: return fs.gas.rho;
        case IOvar::e: return fs.gas.e;
        case IOvar::a: return fs.gas.a;
        case IOvar::velx: return fs.vel.x;
        case IOvar::vely: return fs.vel.y;
        case IOvar::velz: return fs.vel.z;
        default:
            throw runtime_error("Invalid selection for IOvar: "+to_string(i));
        }
        // So we never return from here.
    }

    __host__ __device__
    number estimate_local_dt(Vector3 inorm, Vector3 jnorm, Vector3 knorm, number cfl)
    {
        // We assume that the cells are (roughly) hexagonal and work with
        // velocities normal to the faces.
        number isignal = iLength/(fabs(fs.vel.dot(inorm))+fs.gas.a);
        number jsignal = jLength/(fabs(fs.vel.dot(jnorm))+fs.gas.a);
        number ksignal = kLength/(fabs(fs.vel.dot(knorm))+fs.gas.a);
        return cfl * fmin(fmin(isignal,jsignal),ksignal);
    }

    __host__ __device__
    void eval_dUdt(ConservedQuantities& dUdt, FVFace faces[])
    // These are the spatial (RHS) terms in the semi-discrete governing equations.
    {
        number vol_inv = 1.0/volume;
        auto& fim = faces[face[Face::iminus]];
        auto& fip = faces[face[Face::iplus]];
        auto& fjm = faces[face[Face::jminus]];
        auto& fjp = faces[face[Face::jplus]];
        auto& fkm = faces[face[Face::kminus]];
        auto& fkp = faces[face[Face::kplus]];
        //
        for (int i=0; i < CQI::n; i++) {
            // Integrate the fluxes across the interfaces that bound the cell.
            number surface_integral = fim.area*fim.F[i] - fip.area*fip.F[i]
                + fjm.area*fjm.F[i] - fjp.area*fjp.F[i]
                + fkm.area*fkm.F[i] - fkp.area*fkp.F[i];
            // Then evaluate the derivatives of conserved quantity.
            // Note that conserved quantities are stored per-unit-volume.
            dUdt[i] = vol_inv*surface_integral;
        }
        return;
    } // end eval_dUdt()

}; // end Cell

#endif

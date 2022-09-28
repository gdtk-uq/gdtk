// cell.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef CELL_INCLUDED
#define CELL_INCLUDED

#include <string>
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

    vector<string> names {"iminus", "iplus", "jminus", "jplus", "kminus", "kplus"};
};

int Face_indx_from_name(string name)
{
    if (name == "iminus") return Face::iminus;
    if (name == "iplus") return Face::iplus;
    if (name == "jminus") return Face::jminus;
    if (name == "jplus") return Face::jplus;
    if (name == "kminus") return Face::kminus;
    if (name == "kplus") return Face::kplus;
    throw new runtime_error("Invalid face name: " + name);
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
    // We will keep connections to the pieces compising the cell
    // as indices into global arrays.
    int vtx[8]{0, 0, 0, 0};
    int face[6]{0, 0, 0, 0, 0, 0};

    string toString() {
        string repr = "Cell(pos=" + pos.toString() +
            ", volume=" + to_string(volume) + ", Q=[";
        // for (int i=0; i < TLevels; i++) {
        //     repr += "[";
        //     for (int j=0; j < CQI::n; j++) { repr += to_string(Q[i][j]) + " "; }
        //     repr += "] ";
        // }
        repr += "] )";
        return repr;
    }

    void iovar_set(int i, number val)
    {
        switch (i) {
        case IOvar::posx: pos.x = val; break;
        case IOvar::posy: pos.y = val; break;
        case IOvar::posz: pos.z = val; break;
        case IOvar::p: fs.gas.p = val; break;
        case IOvar::T: fs.gas.T = val; break;
        case IOvar::rho: fs.gas.rho = val; break;
        case IOvar::e: fs.gas.e = val; break;
        case IOvar::a: fs.gas.a = val; break;
        case IOvar::velx: fs.vel.x = val; break;
        case IOvar::vely: fs.vel.y = val; break;
        case IOvar::velz: fs.vel.z = val; break;
        }
    }

    number iovar_get(int i)
    {
        switch (i) {
        case IOvar::posx: return pos.x;
        case IOvar::posy: return pos.y;
        case IOvar::posz: return pos.z;
        case IOvar::p: return fs.gas.p;
        case IOvar::T: return fs.gas.T;
        case IOvar::rho: return fs.gas.rho;
        case IOvar::e: return fs.gas.e;
        case IOvar::a: return fs.gas.a;
        case IOvar::velx: return fs.vel.x;
        case IOvar::vely: return fs.vel.y;
        case IOvar::velz: return fs.vel.z;
        }
    }

    // [TODO] PJ 2022-09-11 construct cell details from the connected faces and vertices.
    // Cell update and IO methods, etc.

}; // end Cell

#endif

// piston.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
//
module piston;

import geom;
import gas;
import gasflow;
import config;
import endcondition;


class Piston {
public:
    size_t indx;
    string label = "";
    double m;   // mass, kg
    double d;   // diameter, m
    double L;   // length, m
    double x;   // position, m
    double vel; // velocity, m/s
    double front_seal_f; // friction factor
    double front_seal_area; // area over which pressure acts
    double back_seal_f;
    double back_seal_area;
    double p_restrain;
    bool is_restrain;
    bool with_brakes;
    bool brakes_on;
    double x_buffer;
    bool hit_buffer;
    EndCondition ecL;
    EndCondition ecR;

    this()
    {
        //
    }
}

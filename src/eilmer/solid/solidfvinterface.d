/**
 * solidfvinterface.d
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-22-04: initial cut
 */

module solidfvinterface;

import geom;

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
    // State
    double T;
    double e;
    double flux;
}




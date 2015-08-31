/**
 * solidfvvertex.d
 *
 * Author: Rowan G. and Peter J.
 * Date: 2015-23-04
 */

module solidfvvertex;

import geom;

class SolidFVVertex {
public:
    size_t id;
    Vector3 pos;
    double areaxy;
    double dTdx;
    double dTdy;

}

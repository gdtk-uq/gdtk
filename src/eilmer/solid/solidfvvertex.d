/**
 * solidfvvertex.d
 *
 * Author: Rowan G. and Peter J.
 * Date: 2015-23-04
 */

module solidfvvertex;

import nm.complex;
import nm.number;
import geom;

class SolidFVVertex {
public:
    size_t id;
    Vector3 pos;
    number areaxy;
    number dTdx;
    number dTdy;
    Vector3[] cloud_pos;
    number*[] cloud_T;

}

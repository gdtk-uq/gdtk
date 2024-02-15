module fsi.damping_models;

import nm.bbla;
import std.math;
// This contains any damping models that people may build for
// the FEM part of the fluid-structure interaction.

void RayleighDampingMatrix(Matrix!double C, Matrix!double K, Matrix!double M, double[] dampingRatio, double[] naturalFreqs) {
    // Build the damping matrix using the Rayleigh 2 parameter, described e.g. in
    // "Hypersonic Fluid-Structure Interaction on a Cantilevered Plate" by
    // Currao et al, equations 4.4 and 5

    // For clarity, unpack the natural frequencies and damping ratios
    double f1 = naturalFreqs[0]; double f2 = naturalFreqs[1];
    double d1 = dampingRatio[0]; double d2 = dampingRatio[1];


    // Eq 5
    double alpha = 2 * f1 * f2 * (d1 * f2 - d2 * f1) / (pow(f2, 2) - pow(f1, 2));
    double beta = 2 * (d2 * f2 - d1 * f1) / (pow(f2, 2) - pow(f1, 2));

    // Eq 4.4
    C._data[] = alpha * M._data[] + beta * K._data[];

    // Modifies C in place, don't return anything
}

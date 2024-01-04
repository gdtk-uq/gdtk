module fsi.kirchofflove;

import std.math;
import std.conv;
import std.algorithm;
import std.stdio;
import core.time;
import std.format;
import std.array;

import nm.number;
import nm.bbla;
import nm.smla;
import fsi;
import geom;

class kirchhoffLovePlate : FEMModel {
public:

    this(string jobName, int id) { super(jobName, id); }

    override void modelSetup() {

        // Set the number of nodes, DoFs and quad points
        nNodes = (myConfig.Nx + 1) * (myConfig.Nz + 1);
        nDoF = nNodes * 3;
        nQuadPoints = myConfig.Nx * myConfig.Nz * 4;

        super.modelSetup();
    }

    override void generateMassStiffnessMatrices() {
        number a = myConfig.length / (2 * myConfig.Nx);
        number b = myConfig.width / (2 * myConfig.Nz);

        Matrix!number KL = LocalStiffnessMatrix(a, b, myConfig.poissonsRatio);
        KL.scale(myConfig.youngsModulus * pow(myConfig.thickness, 3) * a * b / (12 * (1 - pow(myConfig.poissonsRatio, 2))));
        Matrix!number ML = LocalMassMatrix(a, b);
        ML.scale(myConfig.density * myConfig.thickness * a * b);

        size_t globalNodeIndx, globalNodeIndx_2, globalRowIndx, globalColIndx, localRowIndx, localColIndx;
        foreach (k; 0 .. myConfig.Nz) {
            foreach (i; 0 .. myConfig.Nx) {
                foreach (node; 0 .. 4) {
                    globalNodeIndx = LocalNodeToGlobalNode(node, i, k);
                    foreach (DoF; 0 .. 3) {
                        localRowIndx = node * 3 + DoF;
                        globalRowIndx = globalNodeIndx * 3 + DoF;
                        foreach (n_node; 0 .. 4) {
                            globalNodeIndx_2 = LocalNodeToGlobalNode(n_node, i, k);
                            foreach (D_DoF; 0 .. 3) {
                                localColIndx = n_node * 3 + D_DoF;
                                globalColIndx = globalNodeIndx_2 * 3 + D_DoF;
                                M[globalRowIndx, globalColIndx] += ML[localRowIndx, localColIndx];
                                K[globalRowIndx, globalColIndx] += KL[localRowIndx, localColIndx];
                            } // end foreach D_DoF
                        } // end foreach n_node
                    } // end foreach DoF
                } // end foreach node
            } // end foreach i
        } // end foreach k

        foreach (zeroedIndx; zeroedIndices) {
            foreach (DoF; 0 .. nDoFs) {
                if (zeroedIndx == DoF) {
                    M[zeroedIndx, DoF] = 1.0;
                    K[zeroedIndx, DoF] = 1.0;
                } else {
                    M[zeroedIndx, DoF] = 0.0;
                    K[zeroedIndx, DoF] = 0.0;
                    M[DoF, zeroedIndx] = 0.0;
                    K[DoF, zeroedIndx] = 0.0;
                } // end if
            } // end foreach DoF
        } // end foreach zeroedIndx
    } // end GenerateMassStiffnessMatrices

    Matrix!number LocalStiffnessMatrix(number a, number b, number v) {
        // Allocate memory for the local stiffness matrix
        Matrix!number K = new Matrix!number(12); K.zeros();

        // Build the constituitive matrix as per Eq 5.12 in 
        // "Structural Analysis with the Finite Element Method: Linear Statics, Vol 2"
        // The coefficient E/(1-v^2) is applied later

        Matrix!number D = new Matrix!number(3); D.zeros();
        D[0, 0] = 1; D[0, 1] = v;
        D[1, 0] = v; D[1, 1] = 1;
        D[2, 2] = (1 - v) / 2;

        // Use two point quadrature to evaluate the integral
        double[4] xi_quad  = [-1, 1, 1, -1]; xi_quad[]  *= (1 / sqrt(3));
        double[4] eta_quad = [-1, -1, 1, 1]; eta_quad[] *= (1 / sqrt(3));

        // Allocate memory for the B matrix, which is evaluated at each quad point
        Matrix!number B = new Matrix!number(3, 12);

        // Iterate through the quadrature points
        foreach (xi, eta; zip(xi_quad, eta_quad)) {
            // This is transcribed from Eq. 4.40 in
            // "R for Finite Element Analysis of Size-Dependent Microscale Structures"
            // They are the second and mixed derivatives of the shape functions at each node
            B[0, 0] = 3 * xi * (1 - eta) / (4 * a * a);
            B[0, 1] = (3 * xi - 1) * (1 - eta) / (4 * a);
            B[0, 2] = 0;
            B[0, 3] = 3 * xi * (-1 + eta) / (4 * a * a);
            B[0, 4] = (3 * xi + 1) * (1 - eta) / (4 * a);
            B[0, 5] = 0;
            B[0, 6] = 3 * xi * (1 + eta) / (4 * a * a);
            B[0, 7] = (3 * xi + 1) * (1 + eta) / (4 * a);
            B[0, 8] = 0;
            B[0, 9] = 3 * xi * (1 + eta) / (4 * a * a);
            B[0,10] = (3 * xi - 1) * (1 + eta) / (4 * a);
            B[0,11] = 0;
            B[1, 0] = 3 * eta * (1 - xi) / (4 * b * b);
            B[1, 1] = 0;
            B[1, 2] = (3 * eta - 1) * (1 + xi) / (4 * b);
            B[1, 3] = 3 * eta * (1 + xi) / (4 * b * b);
            B[1, 4] = 0;
            B[1, 5] = (3 * eta - 1) * (1 + xi) / (4 * b);
            B[1, 6] = 3 * eta * (-1 + xi) / (4 * b * b);
            B[1, 7] = 0;
            B[1, 8] = (3 * eta + 1) * (1 + xi) / (4 * b);
            B[1, 9] = 3 * eta * (-1 + xi) / (4 * b * b);
            B[1,10] = 0;
            B[1,11] = (3 * eta + 1) * (1 - xi) / (4 * b);
            B[2, 0] = (4 - 3 * xi * xi - 3 * eta * eta) / (4 * a * b);
            B[2, 1] = (-3 * xi * xi + 2 * xi + 1) / (4 * b);
            B[2, 2] = (-3 * eta * eta + 2 * eta + 1) / (4 * a);
            B[2, 3] = (-4 + 3 * xi * xi + 3 * eta * eta) / (4 * a * b);
            B[2, 4] = (-3 * xi * xi - 2 * xi + 1) / (4 * b);
            B[2, 5] = (3 * eta * eta - 2 * eta - 1) / (4 * a);
            B[2, 6] = (4 - 3 * xi * xi - 3 * eta * eta) / (4 * a * b);
            B[2, 7] = (3 * xi * xi + 2 * xi - 1) / (4 * b);
            B[2, 8] = (3 * eta * eta + 2 * eta - 1) / (4 * a);
            B[2, 9] = (-4 + 3 * xi * xi + 3 * eta * eta) / (4 * a * b);
            B[2,10] = (3 * xi * xi - 2 * xi - 1) / (4 * b);
            B[2,11] = (-3 * eta * eta - 2 * eta + 1) / (4 * a);

            Matrix!number BDB = dot(dot(transpose(B), C), B);
            KL._data[] += BDB._data[];
        } // end foreach xi, eta
        return KL;
    } // LocalStiffnessMatrix

    Matrix!number LocalMassMatrix(number a, number b) {
        // Allocate memory for the local mass matrix
        Matrix!number ML = new Matrix!number(12); ML.zeros();

        // Here, we only need the shape functions and their first derivatives,
        // so it's easy enough to write out by hand. So we will loop through the nodes,
        // before looping through the quadrature points to form the N matrices
        // Node locations, in the order specified in Fig. 5.8 of
        // "Structural Analysis with the Finite Element Method: Linear Statics, Vol 2"
        double[4] xn = [-1, 1, 1, -1]; double[4] yn = [-1, -1, 1, 1];

        // Quadrature locations
        double[4] xi_quad = xn[] * (1 / sqrt(3)); double[4] eta_quad = yn[] * (1 / sqrt(3));

        // Allocate memory for the N matrix
        Matrix!number N = new Matrix!number(3, 12);

        // Iterate through quadrature points
        foreach (xi, eta; zip(xi_quad, eta_quad)) {
            // Iterate through nodes on the element
            foreach (i; 0 .. 4) {
                double x = xn[i]; double y = yn[i];
                N[0, i * 3] = (1 + x * xi) * (1 + y * eta) * (2 + x * xi + y * eta - xi * xi - eta * eta) / 8;
                N[0, i * 3 + 1] = a * (xi * xi - 1) * (xi + x) * (1 + y * eta) / 8;
                N[0, i * 3 + 2] = b * (eta * eta - 1) * (eta + y) * (1 + x * xi) / 8;
                N[1, i * 3] = (1 + y * eta) * (2 * xi * x * x + x * (y * eta - 3 * xi * xi - eta * eta + 3) - 2 * xi) / 8;
                N[1, i * 3 + 1] = a * (2 * x * xi + 3 * xi * xi - 1) * (1 + y * eta) / 8;
                N[1, i * 3 + 2] = b * x * (eta * eta - 1) * (eta + y) / 8;
                N[2, i * 3] = (1 + x * xi) * (2 * eta * y * y + y * (x * xi - xi * xi - 3 * eta * eta + 3) - 2 * eta) / 8;
                N[2, i * 3 + 1] = a * y * (xi * xi - 1) * (xi + x) / 8;
                N[2, i * 3 + 2] = b * (2 * y * eta + 3 * eta * eta - 1) * (1 + x * xi) / 8;
            } // end foreach i
            ML.add(dot(transpose(N), N));
        } // end foreach xi, eta

        return ML;
    } // end LocalMassMatrix

    Matrix!number updateForceVector() {
        number a = myConfig.length / (2 * myConfig.Nx);
        number b = myConfig.width / (2 * myConfig.Nz);




        

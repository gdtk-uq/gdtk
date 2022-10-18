// gas_test.cu
// Test program for GasModel and GasState functions.
// PJ 2022-10-18
// Build with
// $ nvcc --expt-relaxed-constexpr -o gas_test gas_test.cu
// Run with
// $ ./gas_test
//
#include <iostream>
#include "number.cu"
#include "gas.cu"

using namespace std;

int main()
{
    cout << "Test for GasModel and GasState." << endl;
    cout << "Air at room conditions." << endl;
    GasState gs;
    gs.p = 101.325e3; // Pascals
    gs.T = 300.0; // degrees K
    gs.update_from_pT();
    cout << "gs=" << gs.toString() << endl;
    //
    number mu, k; gs.trans_coeffs(mu, k);
    cout << "mu=" << mu << "Pa.s  k=" << k << "W/mK" <<endl;
    //
    return 0;
}

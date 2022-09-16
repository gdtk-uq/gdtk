// main.cu
// Toy compressible-flow CFD program to try out CUDA programming.
//
// Initial hack adapts the vector addition example from the CUDA workshop
// to look a bit closer to our Puffin CFD code.
//
// PJ 2022-09-09

#include <iostream>
#include <string>

#include "config.cu"
#include "simulate.cu"

using namespace std;


int main(void)
{
    printf("Chicken compressible-flow CFD\n");
    #ifdef CUDA
    cout << "CUDA flavour of program" << endl;
    #else
    cout << "CPU flavour of program." << endl;
    #endif
    //
    Config::job = "job";
    initialize_simulation();
    march_in_time();
    finalize_simulation();
    //
    return 0;
}

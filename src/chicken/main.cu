// main.cu
// Compressible-flow CFD program to try out CUDA programming.
//

#include <iostream>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <chrono>

#include "include/argagg/argagg.hpp"
#include "config.cu"
#include "simulate.cu"

using namespace std;


int main(int argc, char* argv[])
{
    cout << "Chicken compressible-flow CFD" << endl;
    cout << "Revision-id: PUT_REVISION_STRING_HERE" << endl;
#ifdef CUDA
    cout << "CUDA flavour of program" << endl;
#else
    cout << "CPU flavour of program." << endl;
#endif
    auto clock_start = chrono::system_clock::now();
    //
    argagg::parser argparser {{
        {"help", {"-h", "--help"}, "shows help message", 0},
        {"job", {"-j", "--job"}, "string job-name", 1},
        {"tindx", {"-t", "--tindx"}, "integer time-index", 1},
        {"verbosity", {"-v", "--verbosity"}, "integer level of verbosity for console messages", 1},
    }};
    argagg::parser_results args;
    try {
        args = argparser.parse(argc, argv);
    } catch (const exception& e) {
        cerr << e.what() << endl;
        return 1;
    }
    if (args["help"]) {
        cerr << "chkn-run --job=<string> [--tindx=<int>] [--verbosity=<int>]" << endl;
        return 0;
    }
    Config::verbosity = 0;
    if (args["verbosity"]) Config::verbosity = args["verbosity"].as<int>(0);
    Config::job = "";
    if (args["job"]) Config::job = args["job"].as<string>("");
    if (Config::verbosity > 0) cout << "job name: " << Config::job << endl;
    if (Config::job.length() == 0) {
        cerr << "You need to provide a job name." << endl;
        return 1;
    }
    int tindx_start = 0;
    if (args["tindx"]) tindx_start = args["tindx"].as<int>(0);
    if (Config::verbosity > 0) cout << "tindx_start: " << tindx_start << endl;
    try {
        initialize_simulation(tindx_start);
#ifdef CUDA
        if (!filesystem::exists(filesystem::status("/proc/driver/nvidia"))) {
            throw runtime_error("Cannot find NVIDIA driver in /proc/driver.");
        }
        cudaGetDeviceCount(&Config::nDevices);
        cout << "CUDA devices found: " << Config::nDevices << endl;
        if (Config::nDevices > 0) { march_in_time_using_gpu(); }
#else
        march_in_time_using_cpu_only();
#endif
        finalize_simulation();
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
    }
    auto clock_now = chrono::system_clock::now();
    auto clock_ms = chrono::duration_cast<chrono::milliseconds>(clock_now - clock_start);
    cout << "Done in " << fixed << setprecision(3) << clock_ms.count()/1000.0 << "s" << endl;
    return 0;
}

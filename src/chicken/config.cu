// config.cu
// Include file for chicken: configuration variables.
// PJ 2022-09-17

#ifndef CONFIG_INCLUDED
#define CONFIG_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <stdexcept>
#include "include/nlohmann/json.hpp"

#include "number.cu"
#include "flow.cu"
#include "cell.cu"
#include "block.cu"

using namespace std;
using json = nlohmann::json;

namespace BCCode {
    // Boundary condition codes, to decide what to do for the ghost cells.
    // Periodic boundary conditions should just work if we wrap the index in each direction.
    // There's not enough information here to have arbitrary block connections.
    constexpr int wall_with_slip = 0;
    constexpr int wall_no_slip = 1;
    constexpr int exchange = 2;
    constexpr int inflow = 3;
    constexpr int outflow = 4;

    array<string,5> names{"wall_with_slip", "wall_no_slip", "exchange", "inflow", "outflow"};
};

int BC_code_from_name(string name)
{
    if (name == "wall_with_slip") return BCCode::wall_with_slip;
    if (name == "wall_no_slip") return BCCode::wall_no_slip;
    if (name == "exchange") return BCCode::exchange;
    if (name == "inflow") return BCCode::inflow;
    if (name == "outflow") return BCCode::outflow;
    return BCCode::wall_with_slip;
}


struct BConfig {
    int i, j, k;
    int bcCodes[6];
    int bc_fs[6];
    bool active;

    string toString() {
        ostringstream repr;
        repr << "BConfig(i=" << i << ", j=" << j << ", k=" << k;
        repr << ", bcCodes=["; for (auto c: bcCodes) repr << c << ","; repr << "]";
        repr << ", bc_fs=["; for (auto f : bc_fs) repr << f << ","; repr << "]";
        repr << ", active=" << active;
        repr << ")";
        return repr.str();
    }
}; // end struct BConfig


struct Schedule {
    vector<number> t_change; // times at which the value changes
    vector<number> values; // the corresponding values

    string toString() {
        ostringstream repr;
        repr << "Schedule(";
        repr << "t_change=["; for (auto t : t_change) repr << t << ","; repr << "]";
        repr << "values=["; for (auto v : values) repr << v << ","; repr << "]";
        repr << ")";
        return repr.str();
    }

    number get_value(number t)
    {
        // Select one of our tabulated schedule of values.
        int i = t_change.size() - 1;
        while ((i > 0) && (t < t_change[i])) { i--; }
        return values[i];
    }

    number interpolate_value(number t)
    {
        // Attempt an interpolation of the tabulated schedule of values.
        if (t <= t_change[0]) { return values[0]; }
        int n = t_change.size();
        if (t >= t_change[n-1]) { return values[n-1]; }
        // If we get to this point, we must have at least 2 values in our schedule
        // and we can interpolate between a pair of them.
        int i = n - 1;
        while ((i > 0) && (t < t_change[i])) { i--; }
        number frac = (t-t_change[i])/(t_change[i+1]-t_change[i]);
        number value = (1.0-frac)*values[i] + frac*values[i+1];
        return value;
    }
}; // end struct Schedule


namespace Config {
    int verbosity = 0;
    string job = "job";
    string title = "";
    vector<FlowState> flow_states;
    int nFluidBlocks;
    int nib;
    int njb;
    int nkb;
    vector<int> nics, njcs, nkcs;
    vector<vector<vector<int> > >blk_ids;
    vector<BConfig> blk_configs;
    //
    int print_count = 20;
    int cfl_count = 10;
    Schedule dt_plot_schedule;
    Schedule cfl_schedule;
    number dt_init = 1.0e-6;
    number max_time = 1.0e-3;
    int max_step = 100.0;
    //
    int x_order = 2;
    int t_order = 2;
}

void read_config_file(string fileName)
{
    if (Config::verbosity > 0) cout << "Reading JSON config file." << endl;
    json jsonData;
    ifstream f(fileName, ifstream::binary);
    if (f) {
        f.seekg(0, f.end);
        int length = f.tellg();
        f.seekg(0, f.beg);
        char* text = new char[length];
        f.read(text, length);
        if (f) {
            jsonData = json::parse(text);
        } else {
            cerr << "Could not read all of config.json file." << endl;
        }
        f.close();
        delete[] text;
        if (Config::verbosity > 0) cout << "Finished reading JSON file." << endl;
    } else {
        throw runtime_error("Could not open ifstream for config.json.");
    }
    Config::title = jsonData["title"].get<string>();
    if (Config::verbosity > 0) cout << "  title: " << Config::title << endl;
    //
    // Check that we have consistent gas model.
    map<string, number> gas_model_data = jsonData["gas_model"].get<map<string, number> >();
    if (!approxEquals(GasModel::g, gas_model_data["gamma"], 1.0e-6) ||
        !approxEquals(GasModel::R, gas_model_data["R"], 1.0e-6) ||
        !approxEquals(GasModel::Cv, gas_model_data["Cv"], 1.0e-6)) {
        cerr << "  gas_model_data: gamma:" << gas_model_data["gamma"]
             << " R:" << gas_model_data["R"]
             << " Cv:" << gas_model_data["Cv"]
             << endl;
        throw runtime_error("Incorrect gas model data");
    }
    //
    // Check that we have consistent names in the flow zip archives.
    vector<string> iovar_names = jsonData["iovar_names"].get<vector<string> >();
    for (int i=0; i < IOvar::n; ++i) {
        if (iovar_names[i] != IOvar::names[i]) {
            cout << "Mismatched iovar: " << iovar_names[i] << IOvar::names[i] << endl;
        }
    }
    //
    int n_flow_states = jsonData["n_flow_states"].get<int>();
    vector<json> flow_states_json = jsonData["flow_states"].get<vector<json> >();
    for (auto fsj: flow_states_json) {
        auto gas = GasState{};
        gas.p=fsj["gas"]["p"].get<number>();
        gas.T=fsj["gas"]["T"].get<number>();
        gas.update_from_pT();
        auto vel = fsj["vel"].get<vector<number> >();
        Config::flow_states.push_back(FlowState{gas, Vector3{vel[0], vel[1], vel[2]}});
    }
    if (Config::verbosity > 0) {
        cout << "  flow_states: [";
        for (auto fs : Config::flow_states) cout << fs.toString() << ",";
        cout << "]" << endl;
    }
    //
    // Block configs appear in the JSON file in order of their definition in the Python input script.
    // This index is the block id and is different to the i,j,k indices into the block array.
    //
    Config::nFluidBlocks = jsonData["n_fluid_blocks"].get<int>();
    Config::nib = jsonData["nib"].get<int>();
    Config::njb = jsonData["njb"].get<int>();
    Config::nkb = jsonData["nkb"].get<int>();
    if (Config::verbosity > 0) {
        cout << "  nib: " << Config::nib << endl;
        cout << "  njb: " << Config::njb << endl;
        cout << "  nkb: " << Config::nkb << endl;
    }
    //
    Config::nics = jsonData["nics"].get<vector<int> >();
    Config::njcs = jsonData["njcs"].get<vector<int> >();
    Config::nkcs = jsonData["nkcs"].get<vector<int> >();
    if (Config::verbosity > 0) {
        cout << "  nics: ["; for (auto i : Config::nics) cout << i << ","; cout << "]" << endl;
        cout << "  njcs: ["; for (auto i : Config::njcs) cout << i << ","; cout << "]" << endl;
        cout << "  nkcs: ["; for (auto i : Config::nkcs) cout << i << ","; cout << "]" << endl;
    }
    //
    Config::blk_ids = jsonData["blk_ids"].get<vector<vector<vector<int> > > >();
    if (Config::verbosity > 0) {
        for (int k=0; k < Config::nkb; k++) {
            for (int j=0; j < Config::njb; j++) {
                for (int i=0; i < Config::nib; i++) {
                    cout << "  blk_ids[" << i << "," << j << "," << k << "]=" << Config::blk_ids[i][j][k] << endl;
                }
            }
        }
    }
    //
    vector<json> fluid_blocks_json = jsonData["fluid_blocks"].get<vector<json> >();
    for (auto blk_json : fluid_blocks_json) {
        BConfig blk_config;
        blk_config.i = blk_json["i"].get<int>();
        blk_config.j = blk_json["j"].get<int>();
        blk_config.k = blk_json["k"].get<int>();
        blk_config.active = blk_json["active"].get<bool>();
        map<string,json> bcs_json = blk_json["bcs"].get<map<string,json> >();
        for (auto name : Face::names) {
            map<string,json> bc = bcs_json[name].get<map<string,json> >();
            int i = Face_indx_from_name(name);
            blk_config.bcCodes[i] = BC_code_from_name(bc["tag"].get<string>());
            blk_config.bc_fs[i] = 0; // Default flowstate index.
            if (blk_config.bcCodes[i] == BCCode::inflow) {
                blk_config.bc_fs[i] = bc["flow_state_index"].get<int>();
            }
        }
        if (Config::verbosity > 0) {
            cout << "  blk=" << Config::blk_configs.size() << " config=" << blk_config.toString() << endl;
        }
        Config::blk_configs.push_back(blk_config);
    }
    if (Config::nFluidBlocks != Config::blk_configs.size()) {
        throw runtime_error("Did not read the correct number of block configurations: "+
                            to_string(Config::nFluidBlocks)+" "+to_string(Config::blk_configs.size()));
    }
    //
    Config::x_order = jsonData["x_order"].get<int>();
    Config::t_order = jsonData["t_order"].get<int>();
    Config::dt_init = jsonData["dt_init"].get<number>();
    Config::max_time = jsonData["max_time"].get<number>();
    Config::max_step = jsonData["max_step"].get<int>();
    vector<number> cfl_times = jsonData["cfl_times"].get<vector<number> >();
    vector<number> cfl_values = jsonData["cfl_values"].get<vector<number> >();
    Config::cfl_schedule = Schedule{cfl_times, cfl_values};
    Config::cfl_count = jsonData["cfl_count"].get<int>();
    Config::print_count = jsonData["print_count"].get<int>();
    vector<number> t_change = jsonData["t_change"].get<vector<number> >();
    vector<number> dt_plot = jsonData["dt_plot"].get<vector<number> >();
    Config::dt_plot_schedule = Schedule{t_change, dt_plot};
    if (Config::verbosity > 0) {
        cout << "  x_order=" << Config::x_order << endl;
        cout << "  t_order=" << Config::t_order << endl;
        cout << "  dt_init=" << Config::dt_init << endl;
        cout << "  max_time=" << Config::max_time << endl;
        cout << "  max_step=" << Config::max_step << endl;
        cout << "  cfl_count=" << Config::cfl_count << endl;
        cout << "  cfl_schedule=" << Config::cfl_schedule.toString() << endl;
        cout << "  print_count=" << Config::print_count << endl;
        cout << "  dt_plot_schedule=" << Config::dt_plot_schedule.toString() << endl;
    }
    return;
} // end read_config_file()

#endif

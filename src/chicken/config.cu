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

struct BConfig {
    int i, j, k;
    int initial_fs;
    int bcCodes[6];
    int bc_fs[6];

    string toString() {
        ostringstream repr;
        repr << "BConfig(i=" << i << ", j=" << j << ", k=" << k;
        repr << ", bcCodes=["; for (auto c: bcCodes) repr << c << ","; repr << "]";
        repr << ", bc_fs=["; for (auto f : bc_fs) repr << f << ","; repr << "]";
        repr << ")";
        return repr.str();
    }
};

namespace Config {
    string job = "job";
    string title = "";
    vector<FlowState> flow_states;
    int nib = 1;
    int njb = 1;
    int nkb = 1;
    vector<int> nics, njcs, nkcs;
    vector<vector<vector<int> > >blk_ids;
    vector<BConfig> blk_configs;
}

void read_config_file(string fileName)
{
    cout << "Reading JSON config file." << endl;
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
    } else {
        throw runtime_error("Could not open ifstream for config.json.");
    }
    Config::title = jsonData["title"].get<string>();
    cout << "  title: " << Config::title << endl;
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
    cout << "  flow_states: [";
    for (auto fs : Config::flow_states) cout << fs.toString() << ",";
    cout << "]" << endl;
    //
    // Block configs appear in the JSON file in order of their definition in the Python input script.
    // This index is the block id and is different to the i,j,k indices into the block array.
    //
    Config::nib = jsonData["nib"].get<int>();
    Config::njb = jsonData["njb"].get<int>();
    Config::nkb = jsonData["nkb"].get<int>();
    cout << "  nib: " << Config::nib << endl;
    cout << "  njb: " << Config::njb << endl;
    cout << "  nkb: " << Config::nkb << endl;
    //
    Config::nics = jsonData["nics"].get<vector<int> >();
    Config::njcs = jsonData["njcs"].get<vector<int> >();
    Config::nkcs = jsonData["nkcs"].get<vector<int> >();
    cout << "  nics: ["; for (auto i : Config::nics) cout << i << ","; cout << "]" << endl;
    cout << "  njcs: ["; for (auto i : Config::njcs) cout << i << ","; cout << "]" << endl;
    cout << "  nkcs: ["; for (auto i : Config::nkcs) cout << i << ","; cout << "]" << endl;
    //
    Config::blk_ids = jsonData["blk_ids"].get<vector<vector<vector<int> > > >();
    for (int k=0; k < Config::nkb; k++) {
        for (int j=0; j < Config::njb; j++) {
            for (int i=0; i < Config::nib; i++) {
                cout << "  blk_ids[" << i << "," << j << "," << k << "]=" << Config::blk_ids[i][j][k] << endl;
            }
        }
    }
    //
    int n_fluid_blocks = jsonData["n_fluid_blocks"].get<int>();
    vector<json> fluid_blocks_json = jsonData["fluid_blocks"].get<vector<json> >();
    for (auto blk_json : fluid_blocks_json) {
        BConfig blk_config;
        blk_config.i = blk_json["i"].get<int>();
        blk_config.j = blk_json["j"].get<int>();
        blk_config.k = blk_json["k"].get<int>();
        blk_config.initial_fs = blk_json["initial_flow_state"].get<int>();
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
        cout << "blk_config: " << blk_config.toString() << endl;
        Config::blk_configs.push_back(blk_config);
    }
    return;
}

#endif

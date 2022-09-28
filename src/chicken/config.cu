// config.cu
// Include file for chicken: configuration variables.
// PJ 2022-09-17

#ifndef CONFIG_INCLUDED
#define CONFIG_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "include/nlohmann/json.hpp"

#include "number.cu"
#include "cell.cu"

using namespace std;
using json = nlohmann::json;

namespace Config {
    string job = "job";
    string title = "";
    int nib = 1;
    int njb = 1;
    int nkb = 1;
    vector<int> nics, njcs, nkcs;
    vector<vector<vector<int> > >blk_ids;
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
        cerr << "Could not open ifstream for config.json." << endl;
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
        throw new runtime_error("Incorrect gas model data");
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
    // [TODO] Turn into a more convenient internal format.
    int n_flow_states = jsonData["n_flow_states"].get<int>();
    vector<json> flow_states = jsonData["flow_states"].get<vector<json> >();
    for (int f=0; f < flow_states.size(); f++) {
        cout << "flow_state json: " << flow_states[f] << endl;
    }
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
    // [TODO] Turn into a more convenient internal format.
    int n_fluid_blocks = jsonData["n_fluid_blocks"].get<int>();
    vector<json> fluid_blocks = jsonData["fluid_blocks"].get<vector<json> >();
    for (int b=0; b < fluid_blocks.size(); b++) {
        cout << "fluid_block json: " << fluid_blocks[b] << endl;
    }
    return;
}

#endif

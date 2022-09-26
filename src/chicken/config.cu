// config.cu
// Include file for chicken: configuration variables.
// PJ 2022-09-17

#ifndef CONFIG_INCLUDED
#define CONFIG_INCLUDED

#include "iostream"
#include "fstream"
#include "string"

#include "number.cu"
#include "include/nlohmann/json.hpp"

using namespace std;
using json = nlohmann::json;

namespace Config {
    string job = "job";
    string title = "";
    int nib = 1;
    int njb = 1;
    int nkb = 1;
    vector<int> nics, njcs, nkcs;
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
    Config::nib = jsonData["nib"].get<int>();
    Config::njb = jsonData["njb"].get<int>();
    Config::nkb = jsonData["nkb"].get<int>();
    Config::nics = jsonData["nics"].get<vector<int> >();
    Config::njcs = jsonData["njcs"].get<vector<int> >();
    Config::nkcs = jsonData["nkcs"].get<vector<int> >();
    //
    cout << "  title: " << Config::title << endl;
    cout << "  nib: " << Config::nib << endl;
    cout << "  njb: " << Config::njb << endl;
    cout << "  nkb: " << Config::nkb << endl;
    cout << "  nics: ["; for (auto i : Config::nics) cout << i << ","; cout << "]" << endl;
    cout << "  njcs: ["; for (auto i : Config::njcs) cout << i << ","; cout << "]" << endl;
    cout << "  nkcs: ["; for (auto i : Config::nkcs) cout << i << ","; cout << "]" << endl;
    //
    return;
}

#endif

/** json_helper.d
 * Some convenience functions to help with parsing JSON values from the config files.
 *
 * Author: Peter J.
 * Initial code: 2015-02-05
 */

module util.json_helper;

import std.json;
import std.file;
import std.conv;

JSONValue readJSONfile(string fileName)
{
    string content;
    try {
        content = readText(fileName);
    } catch (Exception e) {
        string msg = text("Failed to read JSON file: ", fileName);
        msg ~= text(" Message is: ", e.msg);
        throw new Exception(msg);
    }
    JSONValue jsonData;
    try {
        jsonData = parseJSON!string(content);
    } catch (Exception e) {
        string msg = text("Failed to parse JSON from file: ", fileName);
        msg ~= text(" Message is: ", e.msg);
        throw new Exception(msg);
    }
    return jsonData;
} // end readJSONfile()


// TODO: lots of repetition here, use templates.

string getJSONstring(JSONValue jsonData, string key, string defaultValue)
{
    string value;
    try {
        value = to!string(jsonData[key].str);
    } catch (Exception e) {
        value = defaultValue;
    }
    return value;
} // end getJSONstring()

int getJSONint(JSONValue jsonData, string key, int defaultValue)
{
    int value;
    try {
        value = to!int(jsonData[key].integer);
    } catch (Exception e) {
        value = defaultValue;
    }
    return value;
} // end getJSONint()

double getJSONdouble(JSONValue jsonData, string key, double defaultValue)
{
    double value;
    try {
        auto json_val = jsonData[key];
        // We wish to accept value like 0.0 or 0
        if (json_val.type() == JSONType.float_) {
            value = json_val.floating;
        } else if (json_val.type() == JSONType.integer) {
            value = to!double(json_val.integer);
        } else {
            value = to!double(json_val.str);
        }
    } catch (Exception e) {
        value = defaultValue;
    }
    return value;
} // end getJSONdouble()

bool getJSONbool(JSONValue jsonData, string key, bool defaultValue)
{
    bool value;
    try {
        value = jsonData[key].type is JSONType.true_;
    } catch (Exception e) {
        value = defaultValue;
    }
    return value;
} // end getJSONbool()

string[] getJSONstringarray(JSONValue jsonData, string key, string[] defaultValue)
{
    string[] value;
    try {
        auto json_values = jsonData[key].array;
        foreach (json_val; json_values) {
            value ~= to!string(json_val.str);
        }
    } catch (Exception e) {
        value = defaultValue;
    }
    return value;
} // end getJSONstringarray()

int[] getJSONintarray(JSONValue jsonData, string key, int[] defaultValue)
{
    int[] value;
    try {
        auto json_values = jsonData[key].array;
        foreach (json_val; json_values) {
            value ~= to!int(json_val.integer);
        }
    } catch (Exception e) {
        value.length = defaultValue.length;
        value[] = defaultValue[];
    }
    return value;
} // end getJSONintarray()

double[] getJSONdoublearray(JSONValue jsonData, string key, double[] defaultValue)
{
    double[] value;
    try {
        auto json_values = jsonData[key].array;
        foreach (json_val; json_values) {
            // We wish to accept value like 0.0 or 0
            if (json_val.type() == JSONType.float_) {
                value ~= json_val.floating;
            } else if (json_val.type() == JSONType.integer) {
                value ~= to!double(json_val.integer);
            } else {
                value ~= to!double(json_val.str);
            }
        }
    } catch (Exception e) {
        value.length = defaultValue.length;
        value[] = defaultValue[];
    }
    return value;
} // end getJSONdoublearray()



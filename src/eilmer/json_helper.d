/** json_helper.d
 * Some convenience functions to help with parsing JSON values from the config files.
 *
 * Author: Peter J.
 * Initial code: 2015-02-05
 */

module json_helper;

import std.json;
import std.conv;

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
	value = to!double(jsonData[key].floating);
    } catch (Exception e) {
	value = defaultValue;
    }
    return value;
} // end getJSONdouble()

bool getJSONbool(JSONValue jsonData, string key, bool defaultValue)
{
    bool value;
    try {
	value = jsonData[key].type is JSON_TYPE.TRUE;
    } catch (Exception e) {
	value = defaultValue;
    }
    return value;
} // end getJSONbool()

int[] getJSONintarray(JSONValue jsonData, string key, int[] defaultValue)
{
    int[] value;
    try {
	auto json_values = jsonData[key].array;
	foreach (json_val; json_values) {
	    value ~= to!int(json_val.integer);
	}
    } catch (Exception e) {
	value = defaultValue;
    }
    return value;
} // end getJSONintarray()

double[] getJSONdoublearray(JSONValue jsonData, string key, double[] defaultValue)
{
    double[] value;
    try {
	auto json_values = jsonData[key].array;
	foreach (json_val; json_values) {
	    value ~= to!double(json_val.floating);
	}
    } catch (Exception e) {
	value = defaultValue;
    }
    return value;
} // end getJSONdoublearray()

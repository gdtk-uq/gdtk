/**
 * A module for the thermal behaviour of a solid (typically as a function of temperature)
 *
 * Author: RJG and KAD
 * Date: 2024-04-27
 *
 * History:
 *   2024-04-29 Added LinearVariation model (RJG)
 */

module lmr.solid.solidthermalmodel;

import std.stdio : File, writefln;
import std.format : format;
import std.typecons : Tuple;
import std.algorithm : min, max;
import std.conv : to;
import std.array : split;
import std.string : strip;
import std.json;
import std.math;

import json_helper;
import util.lua;
import util.lua_service;
import nm.number;
import ntypes.complex;

import lmr.solid.solidstate : SolidState;

import lmrexceptions;
import lmr.lmrerrors;


interface SolidThermalModel
{
public:
    @nogc
    final void updateEnergy(ref SolidState ss) const
    {
        updateProperties(ss);
        _updateEnergy(ss);
    }
    
    @nogc
    final void updateTemperature(ref SolidState ss) const
    {
        _updateTemperature(ss);
        updateProperties(ss);
    }

    @nogc
    void updateProperties(ref SolidState ss) const;

protected:
    // Subclasses need to implement these.
    @nogc
    void _updateEnergy(ref SolidState ss) const;

    @nogc
    void _updateTemperature(ref SolidState ss) const;
}

class ConstantPropertiesSTM : SolidThermalModel
{
public:

    this(double rho, double k, double Cp)
    {
        m_rho = rho;
        m_k = k;
        m_Cp = Cp;
    }

    this(JSONValue jsonData)
    {
        m_rho = getJSONdouble(jsonData, "rho", -1.0);
        m_k = getJSONdouble(jsonData, "k", -1.0);
        m_Cp = getJSONdouble(jsonData, "Cp", -1.0);
    }

    this(lua_State* L, int tblIdx)
    {
        m_rho = getDouble(L, tblIdx, "rho");
        m_k = getDouble(L, tblIdx, "k");
        m_Cp = getDouble(L, tblIdx, "Cp");
    }

    @nogc
    override void updateProperties(ref SolidState ss) const
    {
        ss.rho = m_rho;
        ss.k = m_k;
        ss.Cp = m_Cp;
    }

protected:
    @nogc
    override void _updateEnergy(ref SolidState ss) const
    {
        ss.e = ss.rho * ss.Cp * ss.T;
    }

    @nogc
    override void _updateTemperature(ref SolidState ss) const
    {
        ss.T = ss.e/(ss.rho * ss.Cp);
    }

private:
    double m_rho;
    double m_k;
    double m_Cp;
}

alias LVModelParams = Tuple!(double, "T", double, "rho", double, "k", double, "Cp");

class LinearVariationSTM : SolidThermalModel
{
public:
    this(double e_ref, LVModelParams min, LVModelParams max)
    {
        m_e_ref = e_ref;
        m_min = min;
        m_max = max;
    }

    this(JSONValue jsonData)
    {
        m_e_ref = getJSONdouble(jsonData, "e_ref", 0.0);
        m_min.T = getJSONdouble(jsonData["min"], "T", -1.0);
        m_min.rho = getJSONdouble(jsonData["min"], "rho", -1.0);
        m_min.k = getJSONdouble(jsonData["min"], "k", -1.0);
        m_min.Cp = getJSONdouble(jsonData["min"], "Cp", -1.0);
        m_max.T = getJSONdouble(jsonData["max"], "T", -1.0);
        m_max.rho = getJSONdouble(jsonData["max"], "rho", -1.0);
        m_max.k = getJSONdouble(jsonData["max"], "k", -1.0);
        m_max.Cp = getJSONdouble(jsonData["max"], "Cp", -1.0);
    }

    this(lua_State* L, int tblIdx)
    {
        m_e_ref = getDouble(L, tblIdx, "e_ref");
        lua_getfield(L, tblIdx, "min");
        m_min.T = getDouble(L, -1, "T");
        m_min.rho = getDouble(L, -1, "rho");
        m_min.k = getDouble(L, -1, "k");
        m_min.Cp = getDouble(L, -1, "Cp");
        lua_pop(L, 1);
        lua_getfield(L, -1, "max");
        m_max.T = getDouble(L, -1, "T");
        m_max.rho = getDouble(L, -1, "rho");
        m_max.k = getDouble(L, -1, "k");
        m_max.Cp = getDouble(L, -1, "Cp");
        lua_pop(L, 1);
    }

    @nogc
    override void updateProperties(ref SolidState ss) const
    {
        // Use a constant value (low or high) at edges of range.
        auto T = min(to!number(m_max.T), max(ss.T, to!number(m_min.T)));
        auto w = (T - m_min.T)/(m_max.T - m_min.T);
        ss.rho = (1.0 - w)*m_min.rho + w*m_max.rho;
        ss.k = (1.0 - w)*m_min.k + w*m_max.k;
        ss.Cp = (1.0 - w)*m_min.Cp + w*m_max.Cp;
    }

protected:
    @nogc
    override void _updateEnergy(ref SolidState ss) const
    {
        auto lambda = (m_max.Cp - m_min.Cp)/(m_max.T - m_min.T);
        auto energy = m_min.rho * ((lambda/2) * (ss.T^^2 - m_min.T^^2) + lambda*m_min.T^^2 - lambda*m_min.T*ss.T + m_min.Cp*(ss.T - m_min.T)) + m_e_ref;
        auto e_min = m_e_ref;
        auto e_max = m_min.rho * ((lambda/2) * (ss.T^^2 - m_min.T^^2) + lambda*m_min.T^^2 - lambda*m_min.T*ss.T + m_min.Cp*(m_max.T - m_min.T)) + m_e_ref;
        if (ss.T < m_min.T) {
            ss.e = e_min + m_min.rho * m_min.Cp * (ss.T - m_min.T);
        } else if (m_min.T <= ss.T && ss.T <= m_max.T) {
            ss.e = energy;
        } else {
            ss.e = e_max + m_min.rho * m_max.Cp * (ss.T - m_max.T);
        }
    }
    @nogc
    override void _updateTemperature(ref SolidState ss) const 
    {
        auto lambda = (m_max.Cp - m_min.Cp)/(m_max.T - m_min.T);
        auto sqrtTrm = sqrt(m_min.Cp^^2 * m_min.rho^^2 + 2 * ss.e * lambda * m_min.rho + 2 * lambda * m_min.rho * m_e_ref);
        auto numer = lambda*m_min.rho*m_min.T - m_min.Cp*m_min.rho + sqrtTrm;
        auto denom = lambda*m_min.rho;
        ss.T = numer/denom;
    }

private:
    double m_e_ref;
    LVModelParams m_min;
    LVModelParams m_max;
}

alias TableEntry = Tuple!(double, "T", double, "k", double, "Cp", double, "e");

class TabulatedPropertiesModelSTM : SolidThermalModel {
public:
    this(double rho, string filename)
    {
        File f;
        try {
            f = File(filename, "r");
        }
        catch (Exception e) {
            string errMsg = "Error while initialising the TabulatedPropertiesModel\n";
            errMsg ~= format("Unable to open file: %s\n", filename);
            lmrErrorExit(LmrError.initialisation, errMsg);
        }
        // Throw away first line header.
        f.readln();
        string line;
        while ((line = f.readln()) != null) {
            auto tokens = line.split(",");
            mTab ~= TableEntry();
            mTab[$-1].T = to!float(strip(tokens[0]));
            mTab[$-1].k = to!float(strip(tokens[1]));
            mTab[$-1].Cp = to!float(strip(tokens[2]));
            mTab[$-1].e = to!float(strip(tokens[3]));
        }
        f.close();

        // Construct a linear variation model for each part of the table
        foreach (i; 1 .. mTab.length) {
            auto min = LVModelParams(mTab[i-1].T, rho, mTab[i-1].k, mTab[i-1].Cp);
            auto max = LVModelParams(mTab[i].T, rho, mTab[i].k, mTab[i].Cp);
            mLvms ~= new LinearVariationSTM(mTab[i-1].e, min, max);
        }
    }

    this(JSONValue jsonData)
    {
        auto rho = getJSONdouble(jsonData, "rho", 0.0);
        auto fname = getJSONstring(jsonData, "filename", "");
        this(rho, fname);     
    }

    this(lua_State *L, int tblIdx)
    {
        auto rho = getDouble(L, tblIdx, "rho");
        auto fname = getString(L, tblIdx, "filename");
        this(rho, fname);
    }

    @nogc
    override void updateProperties(ref SolidState ss) const
    {
        if (isNaN(ss.T.re)) {
            _updateTemperature(ss);
        }
        auto idx = findTemperatureInTable(ss.T.re);
        mLvms[idx].updateProperties(ss);
    }

protected:
    /**
     * Return lower index of table entry.
     *
     * For example:
     *
     *  0:  500.0
     *  1: 1000.0
     *  2: 1500.0
     *
     * For T = 1235.0, returns 1. 
     */
    @nogc
    size_t findTemperatureInTable(double T) const
    {
        if (T > mTab[$-1].T) return mLvms.length-1;

        foreach (i; 1 .. mTab.length) {
            if (T < mTab[i].T) return i-1;
        }

        // We shouldn't get here, but we return 0 if we do.
        return 0;
    }

    @nogc
    size_t findEnergyInTable(double e) const
    {
        if (e > mTab[$-1].e) return mLvms.length - 1;

        foreach (i; 1 .. mTab.length) {
            if (e < mTab[i].e) return i-1;
        }

        return 0;
    }

    @nogc
    override void _updateEnergy(ref SolidState ss) const
    {
        auto idx = findTemperatureInTable(ss.T.re);
        mLvms[idx]._updateEnergy(ss);
    }

    @nogc
    override void _updateTemperature(ref SolidState ss) const
    {
        auto idx = findEnergyInTable(ss.e.re);
        mLvms[idx]._updateTemperature(ss);
    }
    

private:
    TableEntry[] mTab;
    LinearVariationSTM[] mLvms;
        
}

SolidThermalModel initSolidThermalModel(JSONValue jsonData)
{
    string modelType = getJSONstring(jsonData, "type", "");

    SolidThermalModel stm;
    try {
        switch (modelType) {
        case "constant_properties":
            stm = new ConstantPropertiesSTM(jsonData);
            break;
        case "linear_variation":
            stm = new LinearVariationSTM(jsonData);
            break;
        case "tabulated_properties":
            stm = new TabulatedPropertiesModelSTM(jsonData);
            break;
        default:
            string errMsg = format("The solid thermal model '%s' is not available.", modelType);
            throw new LmrInitialisationException(errMsg);
        }
    }
    catch (LmrInitialisationException e) {
        string msg = "Error while trying to initialise a solid thermal model.";
        msg ~= e.msg;
        lmrErrorExit(LmrError.initialisation, msg);
    }

    return stm;    
}

SolidThermalModel initSolidThermalModel(lua_State* L, int tblIdx)
{
    string modelType = getString(L, tblIdx, "type");

    SolidThermalModel stm;
    try {
        switch (modelType) {
        case "constant_properties":
            stm = new ConstantPropertiesSTM(L, tblIdx);
            break;
        case "linear_variation":
            stm = new LinearVariationSTM(L, tblIdx);
            break;
        case "tabulated_properties":
            stm = new TabulatedPropertiesModelSTM(L, tblIdx);
            break;
        default:
            string errMsg = format("The solid thermal model '%s' is not available.", modelType);
            throw new LmrInitialisationException(errMsg);
        }
    }
    catch (LmrPreProcessingException e) {
        string msg = "Error while trying to initialise a solid thermal model.";
        msg ~= e.msg;
        lmrErrorExit(LmrError.preprocessing, msg);
    }

    return stm;    
}



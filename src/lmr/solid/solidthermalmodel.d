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

import std.format : format;
import std.typecons : Tuple;
import std.algorithm : min, max;
import std.conv : to;
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



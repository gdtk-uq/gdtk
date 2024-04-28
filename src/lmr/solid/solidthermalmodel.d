/**
 * A module for the thermal behaviour of a solid (typically as a function of temperature)
 *
 * Author: RJG and KAD
 * Date: 2024-04-27
 */

module lmr.solid.solidthermalmodel;

import std.format : format;
import std.json;
import json_helper;
import util.lua;
import util.lua_service;

import lmr.solid.solidstate : SolidState;

import lmrexceptions;
import lmr.lmrerrors;


interface SolidThermalModel {
public:
    @nogc
    final void updateEnergy(ref SolidState ss) const
    {
        updateProperties(ss);
        ss.e = ss.rho * ss.Cp * ss.T;
    }
    
    @nogc
    final void updateTemperature(ref SolidState ss) const
    {
        ss.T = ss.e/(ss.rho * ss.Cp);
        updateProperties(ss);
    }

    @nogc
    void updateProperties(ref SolidState ss) const;
}

class ConstantPropertiesSTM : SolidThermalModel {
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

private:
    double m_rho;
    double m_k;
    double m_Cp;
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



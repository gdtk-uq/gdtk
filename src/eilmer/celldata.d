module celldata;

import std.algorithm.mutation;
import std.math;
import std.conv;
import std.math;
import std.algorithm.mutation : fill;
import std.stdio;

import gas;
import globalconfig;
import fvcell;
import geom;
import fvcell;
import simcore;
import globalconfig;
import nm.complex;
import fluidblockio_old;


class VariableAccess
// VariableAccess gives access to members of FVCell in a way that allows list 
// based access patterns. This means we can generate a list of VariableAccess
// objects and then iterate over them to access FVCell data in a specified 
// order.
{
    public:

    size_t idx;
    size_t num_return = 1;
    double[] buffer;
    

    this(){buffer.length = num_return;}

    @nogc double[] get(FVCell cell){ throw new Error("how did we get here?");}
    @nogc void set(FVCell cell, const double[] value){throw new Error("how did we get here?");}

    @nogc abstract string description();

}

class AuxCellData
// If some data could be added to FVCell, but it isn't a core flow value, it 
// belongs in one of these. This class provides indexing support such that 
// we always know where to find data of a specific type.
{
    public:
    size_t index;

    this(){}
    
    abstract void init(LocalConfig myConfig);
    abstract @nogc void update(FVCell cell, double dt, double time, size_t step);

    static int get_order(string name = "")
    // calculate the order in which these objects will be added
    {
        int counter = 0;
        if (name == CellData.tag) return counter;

        if (GlobalConfig.do_flow_average) {
            counter += 1;
            if (name == FlowAverage.tag) return counter;
        }

        if (GlobalConfig.viscous && GlobalConfig.save_viscous_gradients) {
            counter += 1;
            if (name == CellGradientData.tag) return counter;
        }

        if (GlobalConfig.do_temporal_DFT) {
            counter += 1;
            if (name == GeneralDFT.tag) return counter;
        }

        if (GlobalConfig.solve_electric_field) {
            counter += 1;
            if (name == FieldData.tag) return counter;
        }

        if (name == "") {
            return counter; // hack to get number of items
        } else {
            return -1;
        }
    }

    static AuxCellData[] get_aux_cell_data_items(LocalConfig config)
    // generate a list of items with order determined by get_order()
    {
        AuxCellData[] aux_data;

        aux_data.length = get_order() + 1;

        int i;

        i = get_order(CellData.tag);
        if (i >= 0) aux_data[i] = new CellData();

        i = get_order(FlowAverage.tag);
        if (i >= 0) aux_data[i] = new FlowAverage();

        i = get_order(CellGradientData.tag);
        if (i >= 0) aux_data[i] = new CellGradientData();

        i = get_order(GeneralDFT.tag);
        if (i >= 0) aux_data[i] = new GeneralDFT();

        i = get_order(FieldData.tag);
        if (i >= 0) aux_data[i] = new FieldData();


        foreach (ref aux; aux_data) {
            aux.init(config);
        }

        return aux_data;
    }

}

//=============================================================================
// CellData
//=============================================================================

// Generate all of the accessors for standard FVCell members

template GenCellVariableAccess(string name, string location)
{
    const char[] GenCellVariableAccess = 
    "class "~name~" : VariableAccess
    {
        public:
        @nogc override final double[] get(FVCell cell) {
            buffer[0] = cell."~location~".re; return buffer;
        }
        @nogc override final void set(FVCell cell, const double[] value) { 
            cell."~location~".re = value[0];
        }
        @nogc override final string description(){
            return \"[cell."~location~".re]\";
        }
    }";
}

template GenCellArrayVariableAccess(string name, string location)
{
    const char[] GenCellArrayVariableAccess = 
    "class "~name~" : VariableAccess
    {
        public:
        this(const size_t idx){
            buffer.length = num_return; 
            this.idx = idx;
        }
        @nogc override final double[] get(FVCell cell) {
            buffer[0] = cell."~location~"[idx].re; 
            return buffer;
        }
        @nogc override final void set(FVCell cell, const double[] value) {  
            cell."~location~"[idx].re = value[0];
        }
        @nogc override final string description(){
            return \"[cell."~location~"[idx].re]\";
        }
        private:
        size_t idx;
    }";
}

mixin(GenCellVariableAccess!("AccessPosX", "pos[0]._p[0]"));
mixin(GenCellVariableAccess!("AccessPosY", "pos[0]._p[1]"));
mixin(GenCellVariableAccess!("AccessPosZ", "pos[0]._p[2]"));
mixin(GenCellVariableAccess!("AccessVolume", "volume[0]"));

mixin(GenCellVariableAccess!("AccessVelX", "fs.vel._p[0]"));
mixin(GenCellVariableAccess!("AccessVelY", "fs.vel._p[1]"));
mixin(GenCellVariableAccess!("AccessVelZ", "fs.vel._p[2]"));
mixin(GenCellVariableAccess!("AccessMuT", "fs.mu_t"));
mixin(GenCellVariableAccess!("AccessKT", "fs.k_t"));
mixin(GenCellVariableAccess!("AccessS", "fs.S"));
version(MHD) {
    mixin(GenCellVariableAccess!("AccessMagX", "fs.B._p[0]"));
    mixin(GenCellVariableAccess!("AccessMagY", "fs.B._p[1]"));
    mixin(GenCellVariableAccess!("AccessMagZ", "fs.B._p[2]"));
    mixin(GenCellVariableAccess!("AccessDivB", "fs.divB"));
    mixin(GenCellVariableAccess!("AccessPsi", "fs.psi"));
}

mixin(GenCellVariableAccess!("AccessRho", "fs.gas.rho"));
mixin(GenCellVariableAccess!("AccessQuality", "fs.gas.quality"));
mixin(GenCellVariableAccess!("AccessP", "fs.gas.p"));
mixin(GenCellVariableAccess!("AccessA", "fs.gas.a"));
mixin(GenCellVariableAccess!("AccessMu", "fs.gas.mu"));
mixin(GenCellVariableAccess!("AccessK", "fs.gas.k"));
mixin(GenCellVariableAccess!("AccessU", "fs.gas.u"));
mixin(GenCellVariableAccess!("AccessT", "fs.gas.T"));

mixin(GenCellVariableAccess!("AccessQRadOrg", "Q_rad_org"));
mixin(GenCellVariableAccess!("AccessFRadOrg", "f_rad_org"));
mixin(GenCellVariableAccess!("AccessQRERad", "Q_rE_rad"));

mixin(GenCellVariableAccess!("AccessDtLocal", "dt_local"));

version(multi_T_gas) {
    mixin(GenCellVariableAccess!("AccessDtTherm", "dt_therm"));
    mixin(GenCellArrayVariableAccess!("AccessKModes", "fs.gas.k_modes"));
    mixin(GenCellArrayVariableAccess!("AccessUModes", "fs.gas.u_modes"));  
    mixin(GenCellArrayVariableAccess!("AccessTModes", "fs.gas.T_modes")); 
}
version(multi_species_gas) {
    mixin(GenCellVariableAccess!("AccessDtChem", "dt_chem"));
    mixin(GenCellArrayVariableAccess!("AccessMassF", "fs.gas.massf"));  
}

version(turbulence) {
    mixin(GenCellArrayVariableAccess!("AccessTurb", "fs.turb"));  
}

class CellData : AuxCellData
// An empty auxiliary data item that acts as a pass-through to allow access
// to all of the legacy flow variables while maintaining the same access 
// patterns as all the other I/O
{
    public:

    static tag = "field";

    this(){
        index = AuxCellData.get_order(tag);
    }

    override void init(LocalConfig myConfig){}

    override @nogc void update(FVCell cell, double dt, double time, size_t step){}

    static VariableAccess[string] get_accessors(LocalConfig myConfig)
    {
        VariableAccess[string] acc;

        acc["pos.x"] = new AccessPosX();
        acc["pos.y"] = new AccessPosY();
        acc["pos.z"] = new AccessPosZ();
        acc["volume"] = new AccessVolume();

        acc["rho"] = new AccessRho();
        acc["vel.x"] = new AccessVelX();
        acc["vel.y"] = new AccessVelY();
        acc["vel.z"] = new AccessVelZ();

        version(MHD) {
            acc["B.x"] = new AccessMagX();
            acc["B.y"] = new AccessMagY();
            acc["B.z"] = new AccessMagZ();
            acc["divB"] = new AccessDivB();
            if (myConfig.MHD && myConfig.divergence_cleaning) {
            acc["psi"] = new AccessPsi();
            }
        }

        if (myConfig.include_quality) { 
        acc["quality"] = new AccessQuality();
        }
        acc["p"] = new AccessP();
        acc["a"] = new AccessA();
        acc["mu"] = new AccessMu();
        acc["k"] = new AccessK();

        version(multi_T_gas) {
            foreach(it; 0 .. myConfig.gmodel.n_modes) {
                acc[k_modesName(it)] = new AccessKModes(it);
            }
        }

        acc["mu_t"] = new AccessMuT();
        acc["k_t"] = new AccessKT();
        acc["S"] = new AccessS();

        if (myConfig.radiation) {
            acc["Q_rad_org"] = new AccessQRadOrg();
            acc["f_rad_org"] = new AccessFRadOrg();
            acc["Q_rE_rad"] = new AccessQRERad();
        }

        version(turbulence) {
            foreach(it; 0 .. myConfig.turb_model.nturb) {
                acc[myConfig.turb_model.primitive_variable_name(it)] = new AccessTurb(it);
            }
        }

        version(multi_species_gas) {
            foreach(it; 0 .. myConfig.gmodel.n_species) {
                acc[massfName(myConfig.gmodel, it)] = new AccessMassF(it);
            }
        }

        version(multi_species_gas) {
        if (myConfig.gmodel.n_species > 1) { 
            acc["dt_chem"] = new AccessDtChem();
        }
        }

        acc["u"] = new AccessU();
        acc["T"] = new AccessT();

        version(multi_T_gas) {
            foreach(it; 0 .. myConfig.gmodel.n_modes) {
                acc[u_modesName(it)] = new AccessUModes(it);
                acc[T_modesName(it)] = new AccessTModes(it);
            }
        }
        
        version(multi_T_gas) {
        if (myConfig.gmodel.n_modes > 0) {
            acc["dt_therm"] = new AccessDtTherm();
        }
        }
        
        if (myConfig.with_local_time_stepping) {
        acc["dt_local"] = new AccessDtLocal();
        }

        return acc;
    }
}


//=============================================================================
// FlowAverage
//=============================================================================

// Generate all of the accessors for flow averages held in an auxiliary object

template GenFlowAverageAccess(string name, string location)
{
    const char[] GenFlowAverageAccess = 
    "class "~name~" : VariableAccess
    {
        public:
        this(const size_t idx){
            buffer.length = num_return; 
            aux_idx = idx;
        }
        @nogc override final double[] get(FVCell cell) {
            FlowAverage fa = cast(FlowAverage) cell.aux_cell_data[aux_idx];
            buffer[0] = fa."~location~".re; 
            return buffer;
        }
        @nogc override final void set(FVCell cell, const double[] value) { 
            FlowAverage fa = cast(FlowAverage) cell.aux_cell_data[aux_idx];
            fa."~location~".re = value[0];
        }
        @nogc override final string description(){
            return \"[cell.aux_cell_data[aux_idx]."~location~".re]\";
        }
        private:
        size_t aux_idx;
    }";
}

template GenFlowAverageArrayAccess(string name, string location)
{
    const char[] GenFlowAverageArrayAccess = 
    "class "~name~" : VariableAccess
    {
        public:
        this(const size_t aux_idx, const size_t arr_idx){
            buffer.length = num_return; 
            this.aux_idx = aux_idx; 
            this.arr_idx=arr_idx;
        }
        @nogc override final double[] get(FVCell cell) {
            FlowAverage fa = cast(FlowAverage) cell.aux_cell_data[aux_idx];
            buffer[0] = fa."~location~"[arr_idx].re; 
            return buffer;
        }
        @nogc override final void set(FVCell cell, const double[] value) {
            FlowAverage fa = cast(FlowAverage) cell.aux_cell_data[aux_idx];
            fa."~location~"[arr_idx].re = value[0];
        }
        @nogc override final string description(){return \"[cell.aux_cell_data[aux_idx]."~location~"[arr_idx].re]\";}
        private:
        size_t aux_idx;
        size_t arr_idx;
    }";
}

mixin(GenFlowAverageAccess!("AverageAccessVelX", "vel._p[0]"));
mixin(GenFlowAverageAccess!("AverageAccessVelY", "vel._p[1]"));
mixin(GenFlowAverageAccess!("AverageAccessVelZ", "vel._p[2]"));
mixin(GenFlowAverageAccess!("AverageAccessMuT", "mu_t"));
mixin(GenFlowAverageAccess!("AverageAccessKT", "k_t"));
mixin(GenFlowAverageAccess!("AverageAccessS", "S"));
version(MHD) {
    mixin(GenFlowAverageAccess!("AverageAccessMagX", "B._p[0]"));
    mixin(GenFlowAverageAccess!("AverageAccessMagY", "B._p[1]"));
    mixin(GenFlowAverageAccess!("AverageAccessMagZ", "B._p[2]"));
    mixin(GenFlowAverageAccess!("AverageAccessDivB", "divB"));
    mixin(GenFlowAverageAccess!("AverageAccessPsi", "psi"));
}

mixin(GenFlowAverageAccess!("AverageAccessRho", "rho"));
mixin(GenFlowAverageAccess!("AverageAccessQuality", "quality"));
mixin(GenFlowAverageAccess!("AverageAccessP", "p"));
mixin(GenFlowAverageAccess!("AverageAccessA", "a"));
mixin(GenFlowAverageAccess!("AverageAccessMu", "mu"));
mixin(GenFlowAverageAccess!("AverageAccessK", "k"));
mixin(GenFlowAverageAccess!("AverageAccessU", "u"));
mixin(GenFlowAverageAccess!("AverageAccessT", "T"));

mixin(GenFlowAverageAccess!("AverageAccessQRadOrg", "Q_rad_org"));
mixin(GenFlowAverageAccess!("AverageAccessFRadOrg", "f_rad_org"));
mixin(GenFlowAverageAccess!("AverageAccessQRERad", "Q_rE_rad"));

mixin(GenFlowAverageAccess!("AverageAccessDtLocal", "dt_local"));

version(multi_T_gas) {
    mixin(GenFlowAverageAccess!("AverageAccessDtTherm", "dt_therm"));
    mixin(GenFlowAverageArrayAccess!("AverageAccessKModes", "k_modes"));
    mixin(GenFlowAverageArrayAccess!("AverageAccessUModes", "u_modes"));  
    mixin(GenFlowAverageArrayAccess!("AverageAccessTModes", "T_modes")); 
}
version(multi_species_gas) {
    mixin(GenFlowAverageAccess!("AverageAccessDtChem", "dt_chem"));
    mixin(GenFlowAverageArrayAccess!("AverageAccessMassF", "massf"));  
}

version(turbulence) {
    mixin(GenFlowAverageArrayAccess!("AverageAccessTurb", "turb"));  
}

class FlowAverage : AuxCellData
// Hold flow average data. For 'true' averages this data will need to be
// post-processed (divide by simulation time)
{
    public:

    static tag = "average";

    double rho;
    Vector3 vel;
    
    version(MHD) {
        Vector3 B;
        double divB;
        double psi;
    }

    double quality;
    double p;
    double a;
    double mu;
    double k;

    version(multi_T_gas) {
        double[] k_modes;
    }

    double mu_t;
    double k_t;
    double S;


    double Q_rad_org;
    double f_rad_org;
    double Q_rE_rad;

    version(turbulence) {
        double[] turb;
    }

    version(multi_species_gas) {
        double[] massf;
        double dt_chem;
    }

    double u;
    double T;

    version(multi_T_gas) {
        double[] u_modes;
        double[] T_modes;
    }
    
    version(multi_T_gas) {
        double dt_therm;
    }
    
    double dt_local;

    this(){
        index = AuxCellData.get_order(tag);
    }

    override void init(LocalConfig myConfig)
    {

        version(multi_T_gas) {
            k_modes.length = myConfig.gmodel.n_modes;
        }

        version(turbulence) {
            turb.length = myConfig.turb_model.nturb;
        }

        version(multi_species_gas) {
            massf.length = myConfig.gmodel.n_species;
        }

        version(multi_T_gas) {
            u_modes.length = myConfig.gmodel.n_modes;
            T_modes.length = myConfig.gmodel.n_modes;
        }


        rho = 0.0;
        vel.clear();
    
        version(MHD) {
            B.clear();
            divB = 0.0;
            psi = 0.0;
        }

        quality = 0.0;
        p = 0.0;
        a = 0.0;
        mu = 0.0;
        k = 0.0;

        version(multi_T_gas) {
            fill(k_modes,0.0);
        }

        mu_t = 0.0;
        k_t = 0.0;
        S = 0.0;


        Q_rad_org = 0.0;
        f_rad_org = 0.0;
        Q_rE_rad = 0.0;

        version(turbulence) {
            fill(turb,0.0);
        }

        version(multi_species_gas) {
            fill(massf,0.0);
            dt_chem = 0.0;
        }

        u = 0.0;
        T = 0.0;

        version(multi_T_gas) {
            fill(u_modes,0.0);
            fill(T_modes,0.0);
        }
        
        version(multi_T_gas) {
            dt_therm = 0.0;
        }
        
        dt_local = 0.0;
    
    }

    override @nogc void update(FVCell cell, double dt, double time, size_t step)
    {

        rho += dt*cell.fs.gas.rho.re;
        vel._p[0] += dt*cell.fs.vel.x;
        vel._p[1] += dt*cell.fs.vel.y;
        vel._p[2] += dt*cell.fs.vel.z;

        version(MHD) {
            B._p[0] += dt*cell.fs.B.x;
            B._p[1] += dt*cell.fs.B.y;
            B._p[2] += dt*cell.fs.B.z;
            divB += dt*cell.fs.divB.re;
            if (cell.myConfig.MHD && cell.myConfig.divergence_cleaning) {
                psi += dt*cell.fs.psi.re;
            }
        }

        if (cell.myConfig.include_quality) { 
            quality += dt*cell.fs.gas.quality.re;
        }

        p += dt*cell.fs.gas.p.re;
        a += dt*cell.fs.gas.a.re;
        mu += dt*cell.fs.gas.mu.re;
        k += dt*cell.fs.gas.k.re;

        version(multi_T_gas) {
            foreach(it; 0 .. cell.myConfig.gmodel.n_modes) {
                k_modes[it] += dt*cell.fs.gas.k_modes[it].re;
            }
        }

        mu_t += dt*cell.fs.mu_t.re;
        k_t += dt*cell.fs.k_t.re;
        S += dt*cell.fs.S.re;

        if (cell.myConfig.radiation) {
            Q_rad_org += dt*cell.Q_rad_org.re;
            f_rad_org += dt*cell.f_rad_org.re;
            Q_rE_rad += dt*cell.Q_rE_rad.re;
        }

        version(turbulence) {
            foreach(it; 0 .. cell.myConfig.turb_model.nturb) {
                turb[it] += dt*cell.fs.turb[it].re;
            }
        }

        version(multi_species_gas) {
            foreach(it; 0 .. cell.myConfig.gmodel.n_species) {
                massf[it] += dt*cell.fs.gas.massf[it].re;
            }
        }

        version(multi_species_gas) {
        if (cell.myConfig.gmodel.n_species > 1) { 

            dt_chem += dt*cell.dt_chem;
        }
        }

        u += dt*cell.fs.gas.u.re;
        T += dt*cell.fs.gas.T.re;

        version(multi_T_gas) {
            foreach(it; 0 .. cell.myConfig.gmodel.n_modes) {
                u_modes[it] += dt*cell.fs.gas.u_modes[it].re;
                T_modes[it] += dt*cell.fs.gas.T_modes[it].re;
            }
        }
        
        version(multi_T_gas) {
        if (cell.myConfig.gmodel.n_modes > 0) {
            dt_therm += dt*cell.dt_therm;
        }
        }
        
        if (cell.myConfig.with_local_time_stepping) {
        dt_local += dt*cell.dt_local;
        }

    }

    static VariableAccess[string] get_accessors(LocalConfig myConfig)
    {
        VariableAccess[string] acc;

        size_t idx = AuxCellData.get_order(tag);

        acc["rho"] = new AverageAccessRho(idx);
        acc["vel.x"] = new AverageAccessVelX(idx);
        acc["vel.y"] = new AverageAccessVelY(idx);
        acc["vel.z"] = new AverageAccessVelZ(idx);

        version(MHD) {
            acc["B.x"] = new AverageAccessMagX(idx);
            acc["B.y"] = new AverageAccessMagY(idx);
            acc["B.z"] = new AverageAccessMagZ(idx);
            acc["divB"] = new AverageAccessDivB(idx);
            if (myConfig.MHD && myConfig.divergence_cleaning) {
                acc["psi"] = new AverageAccessPsi(idx);
            }
        }

        if (myConfig.include_quality) { 
            acc["quality"] = new AverageAccessQuality(idx);
        }
        acc["p"] = new AverageAccessP(idx);
        acc["a"] = new AverageAccessA(idx);
        acc["mu"] = new AverageAccessMu(idx);
        acc["k"] = new AverageAccessK(idx);

        version(multi_T_gas) {
            foreach(it; 0 .. myConfig.gmodel.n_modes) {
                acc[k_modesName(it)] = new AverageAccessKModes(idx,it);
            }
        }

        acc["mu_t"] = new AverageAccessMuT(idx);
        acc["k_t"] = new AverageAccessKT(idx);
        acc["S"] = new AverageAccessS(idx);

        if (myConfig.radiation) {
            acc["Q_rad_org"] = new AverageAccessQRadOrg(idx);
            acc["f_rad_org"] = new AverageAccessFRadOrg(idx);
            acc["Q_rE_rad"] = new AverageAccessQRERad(idx);
        }

        version(turbulence) {
            foreach(it; 0 .. myConfig.turb_model.nturb) {
                acc[myConfig.turb_model.primitive_variable_name(it)] = new AverageAccessTurb(idx,it);
            }
        }

        version(multi_species_gas) {
            foreach(it; 0 .. myConfig.gmodel.n_species) {
                acc[massfName(myConfig.gmodel, it)] = new AverageAccessMassF(idx,it);
            }
        }

        version(multi_species_gas) {
        if (myConfig.gmodel.n_species > 1) { 
            acc["dt_chem"] = new AverageAccessDtChem(idx);
        }
        }

        acc["u"] = new AverageAccessU(idx);
        acc["T"] = new AverageAccessT(idx);

        version(multi_T_gas) {
            foreach(it; 0 .. myConfig.gmodel.n_modes) {
                acc[u_modesName(it)] = new AverageAccessUModes(idx,it);
                acc[T_modesName(it)] = new AverageAccessTModes(idx,it);
            }
        }
        
        version(multi_T_gas) {
        if (myConfig.gmodel.n_modes > 0) {
            acc["dt_therm"] = new AverageAccessDtTherm(idx);
        }
        }
        
        if (myConfig.with_local_time_stepping) {
        acc["dt_local"] = new AverageAccessDtLocal(idx);
        }

        return acc;
    }
}

//=============================================================================
// Flow Gradients
//=============================================================================

// we need a different template for array access into FlowGradients
// as the order of indices is [id][x,y,z] instead of [x,y,z][id]
template GenGradArrayAccess(string name, string location, string d)
{
    const char[] GenGradArrayAccess = 
    "class "~name~" : VariableAccess
    {
        public:
        this(const size_t idx){
            buffer.length = num_return; 
            this.idx = idx;
        }
        @nogc override final double[] get(FVCell cell) {
            buffer[0] = cell."~location~"[idx]["~d~"].re; 
            return buffer;
        }
        @nogc override final void set(FVCell cell, const double[] value) {  
            cell."~location~"[idx]["~d~"].re = value[0];
        }
        @nogc override final string description(){
            return \"[cell."~location~"[idx]["~d~"].re]\";
        }
        private:
        size_t idx;
    }";
}

mixin(GenCellVariableAccess!("AccessDUDX", "grad.vel[0][0]"));
mixin(GenCellVariableAccess!("AccessDUDY", "grad.vel[0][1]"));
mixin(GenCellVariableAccess!("AccessDUDZ", "grad.vel[0][2]"));
mixin(GenCellVariableAccess!("AccessDVDX", "grad.vel[1][0]"));
mixin(GenCellVariableAccess!("AccessDVDY", "grad.vel[1][1]"));
mixin(GenCellVariableAccess!("AccessDVDZ", "grad.vel[1][2]"));
mixin(GenCellVariableAccess!("AccessDWDX", "grad.vel[2][0]"));
mixin(GenCellVariableAccess!("AccessDWDY", "grad.vel[2][1]"));
mixin(GenCellVariableAccess!("AccessDWDZ", "grad.vel[2][2]"));

version(multi_species_gas) {
mixin(GenGradArrayAccess!("AccessDMassFDX", "grad.massf", "0"));
mixin(GenGradArrayAccess!("AccessDMassFDY", "grad.massf", "1"));
mixin(GenGradArrayAccess!("AccessDMassFDZ", "grad.massf", "2"));
}

mixin(GenCellVariableAccess!("AccessDTDX", "grad.T[0]"));
mixin(GenCellVariableAccess!("AccessDTDY", "grad.T[1]"));
mixin(GenCellVariableAccess!("AccessDTDZ", "grad.T[2]"));

version(multi_T_gas) {
mixin(GenGradArrayAccess!("AccessDTModesDX", "grad.T_modes", "0"));
mixin(GenGradArrayAccess!("AccessDTModesDY", "grad.T_modes", "1"));
mixin(GenGradArrayAccess!("AccessDTModesDZ", "grad.T_modes", "2"));
}

version(turbulence) {
mixin(GenGradArrayAccess!("AccessDTurbDX", "grad.turb", "0"));
mixin(GenGradArrayAccess!("AccessDTurbDY", "grad.turb", "1"));
mixin(GenGradArrayAccess!("AccessDTurbDZ", "grad.turb", "2"));
}

class CellGradientData : AuxCellData
// An empty auxiliary data item that acts as a pass-through for accessing
// the viscous flow gradients
{
    public:

    static tag = "gradient";

    this(){
        index = AuxCellData.get_order(tag);
    }

    override void init(LocalConfig myConfig){}

    override @nogc void update(FVCell cell, double dt, double time, size_t step){}

    static VariableAccess[string] get_accessors(LocalConfig myConfig)
    {
        VariableAccess[string] acc;

        acc["du_dx"] = new AccessDUDX();
        acc["du_dy"] = new AccessDUDY();
        acc["dv_dx"] = new AccessDVDX();
        acc["dv_dy"] = new AccessDVDY();

        acc["dT_dx"] = new AccessDTDX();
        acc["dT_dy"] = new AccessDTDY();

        if (myConfig.dimensions == 3) {
            acc["du_dz"] = new AccessDUDZ();
            acc["dv_dz"] = new AccessDVDZ();
            acc["dw_dx"] = new AccessDWDX();
            acc["dw_dy"] = new AccessDWDY();
            acc["dw_dz"] = new AccessDWDZ();
            acc["dT_dz"] = new AccessDTDZ();
        }

        version(multi_species_gas) {
        foreach(it; 0 .. myConfig.gmodel.n_species) {
            acc["d"~massfName(myConfig.gmodel, it)~"_dx"] = new AccessDMassFDX(it);
            acc["d"~massfName(myConfig.gmodel, it)~"_dy"] = new AccessDMassFDY(it);
            if (myConfig.dimensions == 3) acc["d"~massfName(myConfig.gmodel, it)~"_dz"] = new AccessDMassFDZ(it);
        }
        }

        version(multi_T_gas) {
        foreach(it; 0 .. myConfig.gmodel.n_modes) {
            acc["d"~k_modesName(it)~"_dx"] = new AccessDTModesDX(it);
            acc["d"~k_modesName(it)~"_dy"] = new AccessDTModesDY(it);
            if (myConfig.dimensions == 3) acc["d"~k_modesName(it)~"_dz"] = new AccessDTModesDZ(it);
        }
        }

        version(turbulence) {
        foreach(it; 0 .. myConfig.turb_model.nturb) {
            acc["d"~myConfig.turb_model.primitive_variable_name(it)~"_dx"] = new AccessDTurbDX(it);
            acc["d"~myConfig.turb_model.primitive_variable_name(it)~"_dy"] = new AccessDTurbDY(it);
            if (myConfig.dimensions == 3) acc["d"~myConfig.turb_model.primitive_variable_name(it)~"_dz"] = new AccessDTurbDZ(it);
        }
        }

        return acc;
    }
}


//=============================================================================
// FlowDFT
//=============================================================================

// accessor for Discrete Fourier Transform (DFT) data

class DFTAccess : VariableAccess
{
    public:

    this(const size_t idx, size_t j, const size_t size)
    {
        aux_idx = idx;
        local_idx = j;
        num_return = 2*size;
        buffer.length = num_return;
    }
    
    @nogc override final double[] get(FVCell cell) 
    {

        GeneralDFT dft = cast(GeneralDFT) cell.aux_cell_data[aux_idx];

        // interleave our data so we have re,im pairs
        foreach (i; 0 .. num_return/2) {
            buffer[2*i + 0] = dft.DFT_local_real[local_idx][i];
            buffer[2*i + 1] = dft.DFT_local_imag[local_idx][i];
        }

        return buffer;
    }

    @nogc override final void set(FVCell cell, const double[] value) 
    { 
        assert(value.length == num_return);

        GeneralDFT dft = cast(GeneralDFT) cell.aux_cell_data[aux_idx];

        foreach (i; 0 .. num_return/2) {
            dft.DFT_local_real[local_idx][i] = value[2*i + 0];
            dft.DFT_local_imag[local_idx][i] = value[2*i + 1];
        }
    }

    @nogc override final string description()
    {
        return "[cell.aux_cell_data[].local_DFT]";
    }

    private size_t aux_idx, local_idx;
}

class GeneralDFT : AuxCellData
// Calculate and hold the DFTs of flow variables
{
    public:

    static tag = "DFT";

    static string[] DFT_names;
    static VariableAccess[] DFT_targets;

    double[][] DFT_local_real;
    double[][] DFT_local_imag;
    size_t DFT_n_modes;

    this(){
        index = AuxCellData.get_order(tag);
    }

    static void set_DFT_targets()
    {
        if (DFT_targets.length > 0) return;

        // the things we want to get the DFT of
        DFT_targets ~= new AccessP(); DFT_names ~= "p";
        DFT_targets ~= new AccessVelX(); DFT_names ~= "vel.x";
        DFT_targets ~= new AccessVelY(); DFT_names ~= "vel.y";
        // DFT_targets ~= new AccessVelZ(); DFT_names ~= "vel.z";
    }

    override void init(LocalConfig myConfig)
    {
        DFT_n_modes = myConfig.DFT_n_modes;

        set_DFT_targets();

        const size_t n_items = DFT_targets.length;

        DFT_local_real.length = n_items;
        DFT_local_imag.length = n_items;

        foreach (i; 0 .. n_items) {
            DFT_local_real[i].length = myConfig.DFT_n_modes;
            DFT_local_imag[i].length = myConfig.DFT_n_modes;
            fill(DFT_local_real[i], 0.0);
            fill(DFT_local_imag[i], 0.0);
        }   

    }

    override @nogc void update(FVCell cell, double dt, double time, size_t step)
    {

        if (step % GlobalConfig.DFT_step_interval != 0) return;

        size_t DFT_step = step / GlobalConfig.DFT_step_interval - 1;
        // If it's the first step, we should set the values rather than incrementing

        foreach (j, acc; DFT_targets) {
            const double val = acc.get(cell)[0]; // assume only getting single values
            foreach (i; 0..DFT_n_modes) {
                DFT_local_real[j][i] += cos(2 * std.math.PI * i * DFT_step / DFT_n_modes) * val;
                DFT_local_imag[j][i] -= sin(2 * std.math.PI * i * DFT_step / DFT_n_modes) * val;
            }
        }

        
    }

    static VariableAccess[string] get_accessors()
    {
        VariableAccess[string] acc;

        size_t idx = AuxCellData.get_order(tag);
        set_DFT_targets();

        foreach (j, name; DFT_names) {
            acc[name] = new DFTAccess(idx, j, GlobalConfig.DFT_n_modes);
        }

        return acc;
    }
}

//=============================================================================
// Cell Limiter Values
//=============================================================================

mixin(GenCellVariableAccess!("AccessVelxPhi", "gradients.velxPhi"));
mixin(GenCellVariableAccess!("AccessVelyPhi", "gradients.velyPhi"));
mixin(GenCellVariableAccess!("AccessVelzPhi", "gradients.velzPhi"));
mixin(GenCellVariableAccess!("AccessRhoPhi", "gradients.rhoPhi"));
mixin(GenCellVariableAccess!("AccessPPhi", "gradients.pPhi"));
mixin(GenCellVariableAccess!("AccessTPhi", "gradients.TPhi"));
mixin(GenCellVariableAccess!("AccessUPhi", "gradients.uPhi"));
version(multi_species_gas) {
    mixin(GenCellArrayVariableAccess!("AccessMassfPhi", "gradients.massfPhi"));
}
version(multi_T_gas) {
    mixin(GenCellArrayVariableAccess!("AccessUModesPhi", "gradients.u_modesPhi"));
    mixin(GenCellArrayVariableAccess!("AccessTModesPhi", "gradients.T_modesPhi"));
}
version(turbulence) {
    mixin(GenCellArrayVariableAccess!("AccessTurbPhi", "gradients.turbPhi"));
}

class CellLimiterData : AuxCellData
// An empty auxiliary data item that acts as a pass-through for accessing
// the unstructured grid limiters
{
    public:

    static tag = "limiter";

    this(){
        index = AuxCellData.get_order(tag);
    }

    override void init(LocalConfig myConfig){}

    override @nogc void update(FVCell cell, double dt, double time, size_t step){}

    static VariableAccess[string] get_accessors(LocalConfig myConfig)
    {
        VariableAccess[string] acc;

        acc["phi_velx"] = new AccessVelxPhi();
        acc["phi_vely"] = new AccessVelyPhi();
        if (myConfig.dimensions == 3) { acc["phi_velz"] = new AccessVelzPhi(); }
        acc["phi_rho"] = new AccessRhoPhi();
        acc["phi_p"] = new AccessPPhi();
        acc["phi_T"] = new AccessTPhi();
        acc["phi_u"] = new AccessUPhi();
        version(multi_species_gas) {
            foreach(sp; 0 .. myConfig.gmodel.n_species) {
                acc["phi_"~massfName(myConfig.gmodel, sp)] = new AccessMassfPhi(sp);
            }
        }
        version(multi_T_gas) {
            foreach(mode; 0 .. myConfig.gmodel.n_modes) {
                acc["phi_"~u_modesName(mode)] = new AccessUModesPhi(mode);
                acc["phi_"~T_modesName(mode)] = new AccessTModesPhi(mode);
            }
        }
        version(turbulence) {
            foreach(it; 0 .. myConfig.turb_model.nturb) {
                acc["phi_"~myConfig.turb_model.primitive_variable_name(it)] = new AccessTurbPhi(it);
            }
        }

        return acc;
    }
}

//=============================================================================
// Field Variables
//=============================================================================

// Generate all of the accessors for Nick's electromagnetic field calculations
mixin(GenCellVariableAccess!("AccessElectricPotential", "electric_potential"));

class FieldData : AuxCellData
// Pipe for accessing field data stored in the cells
{
    public:

    static tag = "efield"; // FIXME: reclaim the "field" tag for myself

    this(){
        index = AuxCellData.get_order(tag);
    }

    override void init(LocalConfig myConfig){}

    override @nogc void update(FVCell cell, double dt, double time, size_t step){}

    static VariableAccess[string] get_accessors(LocalConfig myConfig)
    {
        VariableAccess[string] acc;

        acc["electric_potential"] = new AccessElectricPotential();

        return acc;
    }
}

// io.d: reading and writing and storage for slf
// @author: NNG

module io;

import std.stdio;
import std.format;
import std.string;
import std.conv;

import dyaml;
import gas;
import ntypes.complex;
import nm.number;
import misc;

struct Config {
    string gas_file_name;
    string reaction_file_name;
    size_t N;
    double D;
    double p;
    double T0;
    double T1;
    double[string] Y0;
    double[string] Y1;
    double targetGRR = 1e-10;
    double lewis_number = 0.75;
}

struct Parameters {
    this(ref Config config, GasModel gm) {
        nsp = gm.n_species;
        neq = nsp+1;

        D = config.D;
        p = config.p;
        N = config.N;
        n = N*neq;
        Z.length = N;
        foreach(i; 1 .. N+1) Z[i-1] = i/(N+1.0);
        dZ = 1.0/(N+1.0);

        chi.length = N;
        foreach(i; 0 .. N) {
            chi[i] = evaluate_scalar_dissipation(D, Z[i]);
        }

        // Boundary Conditions
        T0 = config.T0; T1 = config.T1;
        Y0.length = nsp; foreach(isp; 0 .. nsp) Y0[isp] = 0.0;
        Y1.length = nsp; foreach(isp; 0 .. nsp) Y1[isp] = 0.0;

        foreach(sp,mf; config.Y0){
            Y0[gm.species_index(sp)] = mf;
        }
        foreach(sp,mf; config.Y1){
            Y1[gm.species_index(sp)] = mf;
        }

        U0.length = neq; foreach(isp; 0 .. nsp) U0[isp] = Y0[isp]; U0[nsp] = T0;
        U1.length = neq; foreach(isp; 0 .. nsp) U1[isp] = Y1[isp]; U1[nsp] = T1;

        lewis_number = config.lewis_number;
    }

    size_t nsp;
    size_t neq;
    size_t N;
    size_t n;

    number D;
    number p;
    double dZ;
    double T0;
    double T1;

    number[] chi;
    double[] Z;
    double[] Y0;
    double[] Y1;
    double[] U0;
    double[] U1;

    double lewis_number;
}

void write_solution_to_file(ref const Parameters pm, number[] U, string filename){

    File outfile = File(filename, "wb");
    size_t[1] ibuff; double[1] dbuff; // buffer arrays

    ibuff[0] = pm.nsp;  outfile.rawWrite(ibuff);
    ibuff[0] = pm.neq;  outfile.rawWrite(ibuff);
    ibuff[0] = pm.N;    outfile.rawWrite(ibuff);
    ibuff[0] = pm.n;    outfile.rawWrite(ibuff);

    dbuff[0] = pm.D.re; outfile.rawWrite(dbuff);
    dbuff[0] = pm.p.re; outfile.rawWrite(dbuff);
    dbuff[0] = pm.dZ;   outfile.rawWrite(dbuff);
    dbuff[0] = pm.T0;   outfile.rawWrite(dbuff);
    dbuff[0] = pm.T1;   outfile.rawWrite(dbuff);

    foreach(i; 0 .. pm.N)  { dbuff[0] = pm.Z[i];  outfile.rawWrite(dbuff); }
    foreach(i; 0 .. pm.nsp){ dbuff[0] = pm.Y0[i]; outfile.rawWrite(dbuff); }
    foreach(i; 0 .. pm.nsp){ dbuff[0] = pm.Y1[i]; outfile.rawWrite(dbuff); }
    foreach(i; 0 .. pm.n)  { dbuff[0] = U[i].re;     outfile.rawWrite(dbuff); }
    outfile.close();
    return;
}

void read_solution_from_file(ref const Parameters pm, number[] U, string filename){

    File infile = File(filename, "rb");
    size_t[1] ibuff; double[1] dbuff; // buffer arrays

    infile.rawRead(ibuff); if (pm.nsp != ibuff[0]) throw new Error(format("Sim nsp %d does not match file %d", pm.nsp, ibuff[0]));
    infile.rawRead(ibuff); if (pm.neq != ibuff[0]) throw new Error(format("Sim neq %d does not match file %d", pm.neq, ibuff[0]));
    infile.rawRead(ibuff); if (pm.N   != ibuff[0]) throw new Error(format("Sim N   %d does not match file %d", pm.N,   ibuff[0]));
    infile.rawRead(ibuff); if (pm.n   != ibuff[0]) throw new Error(format("Sim n   %d does not match file %d", pm.n,   ibuff[0]));

    // Let's not bother checking these. It's kind of whatever.
    infile.rawRead(dbuff); //if (pm.D   != dbuff[0]) throw new Error(format("Sim D   %e does not match file %e", pm.D,   dbuff[0]));
    infile.rawRead(dbuff); //if (pm.p   != dbuff[0]) throw new Error(format("Sim p   %e does not match file %e", pm.p,   dbuff[0]));
    infile.rawRead(dbuff); //if (pm.dZ  != dbuff[0]) throw new Error(format("Sim dZ  %e does not match file %e", pm.dZ,  dbuff[0]));
    infile.rawRead(dbuff); //if (pm.T0  != dbuff[0]) throw new Error(format("Sim T0  %e does not match file %e", pm.T0,  dbuff[0]));
    infile.rawRead(dbuff); //if (pm.T1  != dbuff[0]) throw new Error(format("Sim T1  %e does not match file %e", pm.T1,  dbuff[0]));

    // At this stage we actually don't want to override the Parameters struct
    foreach(i; 0 .. pm.N)  { infile.rawRead(dbuff); /* pm.Z[i]  = dbuff[0]; */}
    foreach(i; 0 .. pm.nsp){ infile.rawRead(dbuff); /* pm.Y0[i] = dbuff[0]; */}
    foreach(i; 0 .. pm.nsp){ infile.rawRead(dbuff); /* pm.Y1[i] = dbuff[0]; */}
    foreach(i; 0 .. pm.n)  { infile.rawRead(dbuff); U[i].re = dbuff[0]; }
    infile.close();
    return;
}

void write_solution_to_file(ref const Parameters pm, double[] U, string filename){

    File outfile = File(filename, "wb");
    size_t[1] ibuff; double[1] dbuff; // buffer arrays

    ibuff[0] = pm.nsp;  outfile.rawWrite(ibuff);
    ibuff[0] = pm.neq;  outfile.rawWrite(ibuff);
    ibuff[0] = pm.N;    outfile.rawWrite(ibuff);
    ibuff[0] = pm.n;    outfile.rawWrite(ibuff);

    dbuff[0] = pm.D.re; outfile.rawWrite(dbuff);
    dbuff[0] = pm.p.re; outfile.rawWrite(dbuff);
    dbuff[0] = pm.dZ;   outfile.rawWrite(dbuff);
    dbuff[0] = pm.T0;   outfile.rawWrite(dbuff);
    dbuff[0] = pm.T1;   outfile.rawWrite(dbuff);

    foreach(i; 0 .. pm.N)  { dbuff[0] = pm.Z[i];  outfile.rawWrite(dbuff); }
    foreach(i; 0 .. pm.nsp){ dbuff[0] = pm.Y0[i]; outfile.rawWrite(dbuff); }
    foreach(i; 0 .. pm.nsp){ dbuff[0] = pm.Y1[i]; outfile.rawWrite(dbuff); }
    foreach(i; 0 .. pm.n)  { dbuff[0] = U[i].re;     outfile.rawWrite(dbuff); }
    outfile.close();
    return;
}


void write_log_to_file(string[] log, string filename){
    File file = File(filename, "w");
    foreach(line; log){
        file.write(line);
        file.write("\n");
    }
    file.close();
}

Config read_config_from_file(string filename){
/*
    We expect the config variables to be stored in a YAML file.

    Inputs:
        filename - The yaml file to be read.
    Outputs:
        - A new parameters structure, returned by copying.
*/
    Node data = dyaml.Loader.fromFile(filename).load();

    Config cfg = Config();

    cfg.gas_file_name      = data["gas_file_name"].as!string;
    cfg.reaction_file_name = data["reaction_file_name"].as!string;

    cfg.N  = to!size_t(data["N"].as!string);
    cfg.D  = to!double(data["D"].as!string);
    cfg.p  = to!double(data["p"].as!string);
    cfg.T0 = to!double(data["T0"].as!string);
    cfg.T1 = to!double(data["T1"].as!string);

    foreach(Node nd; data["Y0"].mappingKeys) {
        string s = nd.as!string;
        double Ys= to!double(data["Y0"][s].as!string);
        cfg.Y0[s] = Ys;
    }
    foreach(Node nd; data["Y1"].mappingKeys) {
        string s = nd.as!string;
        double Ys= to!double(data["Y1"][s].as!string);
        cfg.Y1[s] = Ys;
    }

    // Optional parameters
    if ("targetGRR" in data) cfg.targetGRR = to!double(data["targetGRR"].as!string);
    if ("lewis_number" in data) cfg.lewis_number = to!double(data["lewis_number"].as!string);

    return cfg;
}

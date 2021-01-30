//univariate_lut_demo.d
// a demo that tests out univariate_lut.d
import std.stdio;
import std.math;
import std.mathspecial;
import std.algorithm;
import std.string;
import std.conv;
import std.datetime.stopwatch;

import nm.univariate_lut;
import gas.gas_model;
import gas.co2gas_sw;


static double F(double rho){
        CO2GasSW gm = new CO2GasSW();
        return gm.get_esat_rho(rho);
}


void main(){
        double rho_min = 0.001;
        double rho_max = 1500;
        int n = 10000;
        uni_lut e_rho_sat_lut = new uni_lut(&F, rho_min, rho_max, n);
        //----DEMO OF READ AND WRITE-------
        string filename = "e_rho_sat_table.dat";
        e_rho_sat_lut.writeLUT(filename);
        e_rho_sat_lut.writeCode();
        uni_lut e_rho_sat_lut2 = new uni_lut(filename);


        double meanError2 = 0;
        double maxError = 0;
        double e;
        double[] es;
        double rho;
        double[] rhos;
        double n_t = 100000;
        for(int i = 0; i != n_t + 1 ; i++){
                rho = i/n_t*(rho_max - rho_min) + rho_min;
                e = e_rho_sat_lut2.quadlookup(rho)/F(rho) - 1;
                es ~= e; rhos ~= rho;
                meanError2 += e*e;
                if (maxError < fabs(e)) writefln("Error is %s at rho= %s", fabs(e), rho);
                maxError = max(fabs(e), maxError);
        }
        File rhoe_error_sat = File("rhoe_error_sat.dat","w");
        rhoe_error_sat.writeln(rhos);
        rhoe_error_sat.writeln(es);
        meanError2 /= sqrt(n_t);
        writefln("Error for %s divisions for %s < rho < %s", n, rho_min, rho_max);
        writefln("MSE: %s, maxError: %s", meanError2, maxError);



        //SOME TIMING TESTS

        int ncycles = 1000;
        StopWatch sw;
        sw.start();
        for(int i = 0; i != ncycles; i++){
                        rho = i/ncycles*(rho_max - rho_min) + rho_min;
                        F(rho);
        }
        sw.stop();
        long t_eval = sw.peek.total!"usecs";
        sw.start();
        for(int i = 0; i != ncycles; i++){
                rho = i/ncycles*(rho_max - rho_min) + rho_min;
                e_rho_sat_lut2.lookup(rho);
        }
        sw.stop();
        long t_lut = sw.peek.total!"usecs" - t_eval;
        writeln("-------------------------------");
        writefln("Time taken for %s evaluations:", ncycles);
        writefln("Direct Eval: %s usecs; LUT: %s usecs;", t_eval, t_lut);
        writefln("Speed Up Factor: %s", t_eval/t_lut);
}

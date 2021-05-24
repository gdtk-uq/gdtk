// tecplot_writer_classic.d
// A simple structured-block tecplot writer.
//
// Author: Peter J. and Rowan G.
// 2021-01-05 Extracted from postprocess.d

module tecplot_writer_classic;

import std.math;
import std.stdio;
import std.conv;
import std.format;
import std.string;
import std.regex;
import std.algorithm;
import std.bitmanip;
import std.stdint;
import std.range;
import nm.complex;
import nm.number;
import gzip;
import fileutil;
import geom;
import gas;
import globalconfig;
import flowsolution;


void write_Tecplot_file(string jobName, string plotDir, FlowSolution soln, int tindx)
{
    ensure_directory_is_present(plotDir);
    auto t = soln.flowBlocks[0].sim_time;
    auto fName = plotDir~"/"~jobName~format("-%.04d", tindx)~".tec";
    auto fp = File(fName, "w");
    fp.writefln("TITLE=\"Job=%s time= %e\"", jobName, t);
    fp.write("VARIABLES= \"X\", \"Y\", \"Z\"");
    int nCtrdVars = 0;
    foreach (var; soln.flowBlocks[0].variableNames) {
        if ( var == "pos.x" || var == "pos.y" || var == "pos.z" ) continue;
        fp.writef(", \"%s\"", var);
        nCtrdVars++;
    }
    fp.write("\n");
    auto ctrdVarsStr = to!string(iota(4,4+nCtrdVars+1));
    foreach (jb; 0 .. soln.nBlocks) {
        auto flow = soln.flowBlocks[jb];
        auto grid = soln.gridBlocks[jb];
        auto nic = flow.nic; auto njc = flow.njc; auto nkc = flow.nkc;
        auto niv = grid.niv; auto njv = grid.njv; auto nkv = grid.nkv;
        fp.writefln("ZONE I=%d J=%d K=%d DATAPACKING=BLOCK", niv, njv, nkv);
        fp.writefln(" SOLUTIONTIME=%e", t);
        fp.writefln(" VARLOCATION=(%s=CELLCENTERED) T=\"fluid-block-%d\"", ctrdVarsStr, jb);
        fp.writefln("# cell-vertex pos.x");
        foreach (k; 0 .. nkv) {
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    fp.writef(" %e", uflowz(grid[i,j,k].x.re));
                }
                fp.write("\n");
            }
        }
        fp.writefln("# cell-vertex pos.y");
        foreach (k; 0 .. nkv) {
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    fp.writef(" %e", uflowz(grid[i,j,k].y.re));
                }
                fp.write("\n");
            }
        }
        fp.writefln("# cell-vertex pos.z");
        foreach (k; 0 .. nkv) {
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    fp.writef(" %e", uflowz(grid[i,j,k].z.re));
                }
                fp.write("\n");
            }
        }
        foreach (var; flow.variableNames) {
            if ( var == "pos.x" || var == "pos.y" || var == "pos.z" ) continue;
            fp.writefln("# cell-centre %s", var);
            foreach (k; 0 .. nkc) {
                foreach (j; 0 .. njc) {
                    foreach (i; 0 .. nic) {
                        fp.writef(" %e", uflowz(flow[var,i,j,k].re));
                    }
                    fp.write("\n");
                }
            }
        }
    } // end for jb in 0..nBlocks
    fp.close();
}

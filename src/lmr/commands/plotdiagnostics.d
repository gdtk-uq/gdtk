/**
 * Module for producing a plot of the residuals and associated diagnostics
 * for a steady-state simulation.
 *
 * Authors: RJG, KAD, NNG
 * Date: 2024-07-09
 */

module lmr.commands.plotdiagnostics;

import std.getopt;
import std.stdio;
import std.string;
import std.format;
import core.stdc.stdlib : system;
import std.process : environment;
import std.file : copy;

import command;

Command plotDiagnosticsCmd;
string cmdName = "plot-diagnostics";

string termPlotCmdFile = "diagnostics-term.gplot";
string rjgPlotCmdFile = "diagnostics-rjg.gplot";
string kadPlotCmdFile = "diagnostics-kad.py";
string nngPlotCmdFile = "diagnostics-nng.py";

static this()
{
    plotDiagnosticsCmd.main = &main_;
    plotDiagnosticsCmd.description =
`Plot diagnostics for a steady-simulation such as residuals and mass balance.
The diagnostics may be plotted in "live" mode to help with monitoring,
or output to file for a snaphot on how convergence has progressed thus far.
`;
    plotDiagnosticsCmd.shortDescription = "Plot diagnostics on convergence for a steady-state simulation.";
    plotDiagnosticsCmd.helpMsg = format(
`lmr %s [options]

%s

This command only makes sense when config.solver_mode = "steady".
It will look for the presence of the file lmrsim/diagnostics/nk-diagnostics.
If no options are provided, the command will write a Gnuplot-type graph
to disk as file: "lmr-ss-diagnostics.pdf".

options ([+] can be repeated):

 -m, --mode=live|file
     sets the mode of operation, options are:
     + live : plot to screen and keep plot alive with a periodic refresh (Ctrl+C to exit)
     + file : write to file based on progress thus far
     examples:
       --mode=live
       -m file
     default: file

 -s, --style=rjg|kad|nng|term
     sets the style of diagnostics display, options are:
     + rjg : a Gnuplot-type graph
     + kad : a matplotlib-type figure on single graph
     + nng : a matplotlib-type figure split on two graphs
     + term : an ASCII rendering of the plot in the terminal
              (if used with -m file no file is written, but a 
              single plot is displayed and returns)
     examples:
       --style=kad
       -s nng
     default: rjg

 -o, --output=filename.pdf
     sets the name for the output file.
     For --style=rjg, the only allowable extensions are: png and pdf.
     These will be used to control the Gnuplot output.
     For --style=kad|nng, the allowable extension are whatever can
     be handled by matplotlib's savefig() function.
     examples:
       -o diag-plot.png
       --output=diagnostics-fig.pdf
     default: lmr-ss-diagnostics.pdf

 -c, --copy-script
     gets a copy of the plotting script placed in the local directory.
     This is useful if one wants to tweak or tune the output aesthetics
     for their purposes.

  -v, --verbose [+]
     Increase verbosity.
`, cmdName, plotDiagnosticsCmd.shortDescription);
 
}

int main_(string[] args)
{
    int verbosity = 0;
    enum Mode { live, file }
    Mode mode = Mode.file;
    enum Style { term, rjg, kad, nng }
    Style style = Style.rjg;
    string outputFile = "lmr-ss-diagnostics.pdf";
    bool copyScript = false;

    try {
        getopt(args,
               config.bundling,
               "v|verbosity", &verbosity,
               "m|mode", &mode,
               "s|style", &style,
               "o|output", &outputFile,
               "c|copy-script", &copyScript
               );
    }
    catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (verbosity > 0) {
        writefln("lmr %s: Begin plotting.", cmdName);
    }

    // Figure out extension on file, or set as PDF.
    string fileExt = "pdf";
    if (mode == Mode.file) {
        auto tokens = outputFile.split(".");
        if (tokens.length > 0) {
            fileExt = tokens[$-1];
        }
    }
    

    auto shareDir = environment.get("DGD") ~ "/share/";
    string plotCmdFile;
    string plotCmdFileLocalUse;
    string cmd;
    string cmdForLocalUse;
    final switch (style) {
    case Style.term:
        plotCmdFileLocalUse = termPlotCmdFile;
        plotCmdFile = shareDir ~ plotCmdFileLocalUse;
        final switch (mode) {
        case Mode.live:
            cmd = "gnuplot -c " ~ plotCmdFile ~ " live";
            cmdForLocalUse = "gnuplot -c " ~ plotCmdFileLocalUse ~ " live";
            break;
        case Mode.file:
            cmd = "gnuplot " ~ plotCmdFile;
            cmdForLocalUse = "gnuplot " ~ plotCmdFileLocalUse;
            break;
        }
        break;
        
    case Style.rjg:
        plotCmdFileLocalUse = rjgPlotCmdFile;
        plotCmdFile = shareDir ~ plotCmdFileLocalUse;
        final switch (mode) {
        case Mode.live:
            cmd = "gnuplot -c " ~ plotCmdFile ~ " live";
            cmdForLocalUse = "gnuplot -c " ~ plotCmdFileLocalUse ~ " live";
            break;
        case Mode.file:
            switch (fileExt) {
            case "pdf":
                cmd = "gnuplot -c " ~ plotCmdFile ~ " pdf " ~ outputFile;
                cmdForLocalUse = "gnuplot -c " ~ plotCmdFileLocalUse ~ " pdf " ~ outputFile;
                break;
            case "png":
                cmd = "gnuplot -c " ~ plotCmdFile ~ " png " ~ outputFile;
                cmdForLocalUse = "gnuplot -c " ~ plotCmdFileLocalUse ~ " png " ~ outputFile;
                break;
            default:
                string errMsg = format("The extension type '%s' on the output file is unrecognised.\n", fileExt);
                errMsg ~= "Allowable options for 'rjg' style are:\n";
                errMsg ~= "   'pdf'\n";
                errMsg ~= "   'png'\n";
                throw new Error(errMsg);
            }
            break;
        }
        break;

    case Style.kad, Style.nng:
        plotCmdFileLocalUse = (style == Style.kad) ? kadPlotCmdFile : nngPlotCmdFile;
        plotCmdFile = shareDir ~ plotCmdFileLocalUse;
        final switch (mode) {
        case Mode.live:
            cmd = "python " ~ plotCmdFile ~ " live";
            cmdForLocalUse = "python " ~ plotCmdFileLocalUse ~ " live";
            break;
        case Mode.file:
            cmd = "python " ~ plotCmdFile ~ " to-file " ~ outputFile;
            cmdForLocalUse = "python " ~ plotCmdFileLocalUse ~ " to-file " ~ outputFile;
        }
        break;
    }
    

    if (copyScript) {
        writefln("Plotting script copy requested. Copying file: %s", plotCmdFile);
        copy(plotCmdFile, plotCmdFileLocalUse);
        writeln("Command to execute is:");
        writefln("> %s", cmdForLocalUse);
    }

    if (verbosity >= 2) {
        writefln("lmr %s: Execute...", cmdName);
        writeln(cmd);
    }

    // Now call command
    auto flag  = system(cmd.toStringz);

    
    if (verbosity > 0) {
        writefln("lmr %s: Done.", cmdName);
    }

    return flag;
    
}

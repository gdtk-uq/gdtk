/** lmrconfig.d
 * Module for configuration of Eilmer program itself.
 *
 * Authors: RJG, PAJ, KAD, NNG
 * Date: 2022-08-06
 */

module lmr.lmrconfig;

import std.conv : to;
import std.format : format;
import std.json;
import std.process : environment;
import std.stdio;

import util.json_helper : readJSONfile;

import lmr.globalconfig;

struct LmrCfg {
    immutable string simDir;
    immutable string jobFile;
    immutable string cfgFile;
    immutable string ctrlFile;
    immutable string progFile;
    immutable string nkCfgFile;
    immutable string nkCmdsFile;
    immutable string blkIdxFmt;
    immutable string cellIdxFmt;
    immutable string snapshotDir;
    immutable string snapshotIdxFmt;
    immutable string historyDir;
    immutable string historyPrefix;
    immutable int initialFieldDir;
    immutable string fluidMetadataFile;
    immutable string solidMetadataFile;
    immutable string limiterMetadataFile;
    immutable string residualMetadataFile;
    immutable string gradientMetadataFile;
    immutable string restartFile;
    immutable string timesFile;
    immutable string referenceResidualsFile;
    immutable string fluidPrefix;
    immutable string solidPrefix;
    immutable string limiterPrefix;
    immutable string residualPrefix;
    immutable string gradientPrefix;
    immutable string loadsDir;
    immutable string loadsPrefix;
    immutable string gridPrefix;
    immutable string gridDir;
    immutable string gridDirectory;
    immutable string gridMetadataFile;
    immutable string savedSgridDir;
    immutable string gzipExt;
    immutable string blkListFile;
    immutable string vtkDir;
    immutable string mpimapFile;
    immutable string mappedCellsFile;
    immutable string transResidFile;
    immutable string dblVarFmt;
}

LmrCfg lmrCfg;
JSONValue lmrJSONCfg;

static this()
{
    // Read base config into lmrJSONCfg
    readLmrConfig();

    // Populate struct with derived config information
    lmrCfg.simDir = lmrJSONCfg["simulation-directory"].str;
    lmrCfg.jobFile = lmrJSONCfg["job-filename"].str;
    lmrCfg.cfgFile = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["config-filename"].str;
    lmrCfg.ctrlFile = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["control-filename"].str;
    lmrCfg.progFile = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["progress-filename"].str;
    lmrCfg.nkCfgFile = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["newton-krylov-config-filename"].str;
    lmrCfg.nkCmdsFile = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["newton-krylov-commands-filename"].str;
    lmrCfg.blkIdxFmt = lmrJSONCfg["block-index-format"].str;
    lmrCfg.cellIdxFmt = lmrJSONCfg["cell-index-format"].str;
    lmrCfg.snapshotDir = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["snapshot-directory"].str;
    lmrCfg.snapshotIdxFmt = lmrJSONCfg["snapshot-index-format"].str;
    lmrCfg.historyDir = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["history-directory"].str;
    lmrCfg.historyPrefix = lmrJSONCfg["history-prefix"].str;
    lmrCfg.initialFieldDir = to!int(lmrJSONCfg["initial-field-directory"].integer);
    lmrCfg.fluidMetadataFile = lmrCfg.snapshotDir ~ "/" ~ lmrJSONCfg["fluid-prefix"].str ~ lmrJSONCfg["metadata-extension"].str;
    lmrCfg.solidMetadataFile = lmrCfg.snapshotDir ~ "/" ~ lmrJSONCfg["solid-prefix"].str ~ lmrJSONCfg["metadata-extension"].str;
    lmrCfg.limiterMetadataFile = lmrCfg.snapshotDir ~ "/" ~ lmrJSONCfg["limiter-prefix"].str ~ lmrJSONCfg["metadata-extension"].str;
    lmrCfg.residualMetadataFile = lmrCfg.snapshotDir ~ "/" ~ lmrJSONCfg["residual-prefix"].str ~ lmrJSONCfg["metadata-extension"].str;
    lmrCfg.gradientMetadataFile = lmrCfg.snapshotDir ~ "/" ~ lmrJSONCfg["gradient-prefix"].str ~ lmrJSONCfg["metadata-extension"].str;
    lmrCfg.restartFile = lmrCfg.snapshotDir ~ "/" ~ lmrJSONCfg["restart-filename"].str;
    lmrCfg.timesFile = lmrCfg.snapshotDir ~ "/" ~ lmrJSONCfg["times-filename"].str;
    lmrCfg.referenceResidualsFile = lmrCfg.snapshotDir ~ "/" ~ lmrJSONCfg["reference-residuals-file"].str;
    lmrCfg.fluidPrefix = lmrJSONCfg["fluid-prefix"].str;
    lmrCfg.solidPrefix = lmrJSONCfg["solid-prefix"].str;
    lmrCfg.limiterPrefix = lmrJSONCfg["limiter-prefix"].str;
    lmrCfg.residualPrefix = lmrJSONCfg["residual-prefix"].str;
    lmrCfg.gradientPrefix = lmrJSONCfg["gradient-prefix"].str;
    lmrCfg.loadsDir = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["loads-directory"].str;
    lmrCfg.loadsPrefix = lmrJSONCfg["loads-prefix"].str;
    lmrCfg.gridPrefix = lmrJSONCfg["grid-prefix"].str;
    lmrCfg.gridDir = lmrJSONCfg["grid-directory"].str;
    lmrCfg.gridDirectory = lmrCfg.simDir ~ "/" ~ lmrCfg.gridDir;
    lmrCfg.gridMetadataFile = lmrCfg.gridDirectory ~ "/" ~ lmrJSONCfg["grid-prefix"].str ~ lmrJSONCfg["metadata-extension"].str;
    lmrCfg.savedSgridDir = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["saved-sgrid-directory"].str;
    lmrCfg.gzipExt = lmrJSONCfg["gzip-extension"].str;
    lmrCfg.blkListFile = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["block-list-filename"].str;
    lmrCfg.vtkDir = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["vtk-output-directory"].str;
    lmrCfg.mpimapFile = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["mpimap-filename"].str;
    lmrCfg.mappedCellsFile = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["mappedcells-filename"].str;
    lmrCfg.transResidFile = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["transient-residuals-filename"].str;
    lmrCfg.dblVarFmt = lmrJSONCfg["double-var-format"].str;


}

/**
 * Read Eilmer program configuration file.
 *
 * Authors: RJG
 * Date: 2022-08-06
 */
void readLmrConfig()
{
    auto lmrCfgFile = environment.get("DGD") ~ "/etc/lmr.cfg";
    lmrJSONCfg = readJSONfile(lmrCfgFile);
}

/**
 * Return the snapshot directory for index 'snapshot'
 *
 * Authors: RJG
 * Date: 2023-06-27
 */
string snapshotDirectory(int snapshot)
{
    return lmrCfg.snapshotDir ~
    	"/" ~
	    format(lmrCfg.snapshotIdxFmt, snapshot);
}


/**
 * Return the fluid solution filename for a single block ('blkId') as a string.
 *
 * Authors: RJG
 * Date: 2023-06-27
 */
string fluidFilename(int snapshot, int blkId)
{
    string fname = lmrCfg.snapshotDir ~
        "/" ~
    	format(lmrCfg.snapshotIdxFmt, snapshot) ~
	    "/" ~
	    lmrCfg.fluidPrefix ~ "-" ~ format(lmrCfg.blkIdxFmt, blkId);
    if (GlobalConfig.field_format == "gziptext")
	    fname ~= lmrCfg.gzipExt;
    return fname;
}

/**
 * Return the solid solution filename for a single block ('blkId') as a string.
 *
 * Authors: RJG
 * Date: 2024-02-25
 */
string solidFilename(int snapshot, int blkId)
{
    string fname = lmrCfg.snapshotDir ~
    	"/" ~
	    format(lmrCfg.snapshotIdxFmt, snapshot) ~
	    "/" ~
	    lmrCfg.solidPrefix ~ "-" ~ format(lmrCfg.blkIdxFmt, blkId);
    if (GlobalConfig.field_format == "gziptext")
	    fname ~= lmrCfg.gzipExt;
    return fname;
}

/**
 * Return the limiter values filename for a single block ('blkId') as a string.
 *
 * Authors: RJG
 * Date: 2023-08-13
 */
string limiterFilename(int snapshot, int blkId)
{
    string fname = lmrCfg.snapshotDir ~
    	"/" ~
	    format(lmrCfg.snapshotIdxFmt, snapshot) ~
	    "/" ~
	    lmrCfg.limiterPrefix ~ "-" ~ format(lmrCfg.blkIdxFmt, blkId);
    if (GlobalConfig.field_format == "gziptext")
	    fname ~= lmrCfg.gzipExt;
    return fname;
}

/**
 * Return the residual values filename for a single block ('blkId') as a string.
 *
 * Authors: KAD and RJG
 * Date: 2024-03-07
 */
string residualFilename(int snapshot, int blkId)
{
    string fname = lmrCfg.snapshotDir ~
        "/" ~
        format(lmrCfg.snapshotIdxFmt, snapshot) ~
        "/" ~
        lmrCfg.residualPrefix ~ "-" ~ format(lmrCfg.blkIdxFmt, blkId);
    if (GlobalConfig.field_format == "gziptext")
        fname ~= lmrCfg.gzipExt;
    return fname;
}

/**
 * Return the gradient values filename for a single block ('blkId') as a string.
 *
 * Authors: KAD and RJG
 * Date: 2024-08-01
 */
string gradientFilename(int snapshot, int blkId)
{
    string fname = lmrCfg.snapshotDir ~
        "/" ~
	    format(lmrCfg.snapshotIdxFmt, snapshot) ~
	    "/" ~
	    lmrCfg.gradientPrefix ~ "-" ~ format(lmrCfg.blkIdxFmt, blkId);
    if (GlobalConfig.field_format == "gziptext")
	    fname ~= lmrCfg.gzipExt;
    return fname;
}

/**
 * Return the loads filename for a single block+boundary combo as a string.
 *
 * Authors: RJG
 * Date: 2023-11-19
 */
string loadsFilename(int tindx, int blkId, int bndryId, string group)
{
    string fname = lmrCfg.loadsDir ~
        "/" ~ format(lmrCfg.snapshotIdxFmt, tindx) ~
        "/" ~
        "blk-" ~ format(lmrCfg.blkIdxFmt, blkId) ~
        "-bndry-" ~ format("%d", bndryId) ~
        "-" ~ group ~
         ".dat";
    return fname;
}
/**
 * Return the grid filename for a single block ('id') as a string.
 *
 * Authors: RJG
 * Date: 2022-08-06
 */
string gridFilename(int snapshot, int blkId)
{
    string gname = lmrCfg.snapshotDir ~
        "/" ~
        format(lmrCfg.snapshotIdxFmt, snapshot) ~
        "/" ~
        lmrCfg.gridPrefix ~ "-" ~ format(lmrCfg.blkIdxFmt, blkId);
    if (GlobalConfig.grid_format == "gziptext") gname ~= lmrCfg.gzipExt;
    return gname;
}

/**
 * Return the name of a history file based on block and cell.
 *
 * Authors: RJG
 * Date: 2024-02-07
 * History:
 *    2024-03-19 -- remove hist prefix,
 *                  we think it's somewhat obvious
 *                  given these files live in a hist/ area
 */
string historyFilename(size_t hcellIdx, size_t blkId, size_t cellId)
{
    string hname = lmrCfg.historyDir ~ "/" ~
         "hc-" ~ format("%02d", hcellIdx) ~
         "-blk-" ~ format(lmrCfg.blkIdxFmt, blkId) ~
         "-cell-" ~ format(lmrCfg.cellIdxFmt, cellId) ~
         ".dat";
    return hname;
}

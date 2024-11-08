/** lmrconfig.d
 * Module for configuration of Eilmer program itself.
 *
 * Authors: RJG, PAJ, KAD, NNG
 * Date: 2022-08-06
 */

module lmrconfig;

import std.stdio;
import std.conv : to;
import std.format : format;
import std.process : environment;
import std.json;
import util.json_helper : readJSONfile;

import globalconfig;

struct LmrCfg {
    shared immutable string simDir;
    shared immutable string jobFile;
    shared immutable string cfgFile;
    shared immutable string ctrlFile;
    shared immutable string progFile;
    shared immutable string nkCfgFile;
    shared immutable string blkIdxFmt;
    shared immutable string cellIdxFmt;
    shared immutable string snapshotDir;
    shared immutable string snapshotIdxFmt;
    shared immutable string historyDir;
    shared immutable string historyPrefix;
    shared immutable int initialFieldDir;
    shared immutable string fluidMetadataFile;
    shared immutable string solidMetadataFile;
    shared immutable string limiterMetadataFile;
    shared immutable string residualMetadataFile;
    shared immutable string gradientMetadataFile;
    shared immutable string restartFile;
    shared immutable string timesFile;
    shared immutable string referenceResidualsFile;
    shared immutable string fluidPrefix;
    shared immutable string solidPrefix;
    shared immutable string limiterPrefix;
    shared immutable string residualPrefix;
    shared immutable string gradientPrefix;
    shared immutable string loadsDir;
    shared immutable string loadsPrefix;
    shared immutable string gridPrefix;
    shared immutable string gridDir;
    shared immutable string gridDirectory;
    shared immutable string gridMetadataFile;
    shared immutable string savedSgridDir;
    shared immutable string gzipExt;
    shared immutable string blkListFile;
    shared immutable string vtkDir;
    shared immutable string mpimapFile;
    shared immutable string mappedCellsFile;
    shared immutable string transResidFile;
    shared immutable string dblVarFmt;
    shared immutable string revisionId = "PUT_REVISION_STRING_HERE";
    shared immutable string fullRevisionId = "PUT_FULL_REVISION_STRING_HERE";
    shared immutable string revisionDate = "PUT_REVISION_DATE_HERE";
    shared immutable string compilerName = "PUT_COMPILER_NAME_HERE";
    shared immutable string buildDate = "PUT_BUILD_DATE_HERE";
};

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

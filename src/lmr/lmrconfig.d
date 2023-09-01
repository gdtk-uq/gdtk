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
import json_helper : readJSONfile;

import globalconfig;

struct LmrCfg {
    immutable string simDir;
    immutable string jobFile;
    immutable string cfgFile;
    immutable string ctrlFile;
    immutable string nkCfgFile;
    immutable string blkIdxFmt;
    immutable string snapshotDir;
    immutable string snapshotIdxFmt;
    immutable int initialFieldDir;
    immutable string flowMetadataFile;
    immutable string limiterMetadataFile;
    immutable string restartFile;
    immutable string referenceResidualsFile;
    immutable string flowPrefix;
    immutable string limiterPrefix;
    immutable string gridPrefix;
    immutable string gridDir;
    immutable string gridDirectory;
    immutable string gzipExt;
    immutable string blkListFile;
    immutable string vtkDir;
    immutable string mpimapFile;
    immutable string dblVarFmt;
    immutable string revisionId = "PUT_REVISION_STRING_HERE";
    immutable string fullRevisionId = "PUT_FULL_REVISION_STRING_HERE";
    immutable string revisionDate = "PUT_REVISION_DATE_HERE";
    immutable string compilerName = "PUT_COMPILER_NAME_HERE";
    immutable string buildDate = "PUT_BUILD_DATE_HERE";
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
    lmrCfg.nkCfgFile = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["newton-krylov-config-filename"].str;
    lmrCfg.blkIdxFmt = lmrJSONCfg["block-index-format"].str;
    lmrCfg.snapshotDir = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["snapshot-directory"].str;
    lmrCfg.snapshotIdxFmt = lmrJSONCfg["snapshot-index-format"].str;
    lmrCfg.initialFieldDir = to!int(lmrJSONCfg["initial-field-directory"].integer);
    lmrCfg.flowMetadataFile = lmrCfg.snapshotDir ~ "/" ~ lmrJSONCfg["flow-prefix"].str ~ lmrJSONCfg["metadata-extension"].str;
    lmrCfg.limiterMetadataFile = lmrCfg.snapshotDir ~ "/" ~ lmrJSONCfg["limiter-prefix"].str ~ lmrJSONCfg["metadata-extension"].str;
    lmrCfg.restartFile = lmrCfg.snapshotDir ~ "/" ~ lmrJSONCfg["restart-file"].str;
    lmrCfg.referenceResidualsFile = lmrCfg.snapshotDir ~ "/" ~ lmrJSONCfg["reference-residuals-file"].str;
    lmrCfg.flowPrefix = lmrJSONCfg["flow-prefix"].str;
    lmrCfg.limiterPrefix = lmrJSONCfg["limiter-prefix"].str;
    lmrCfg.gridPrefix = lmrJSONCfg["grid-prefix"].str;
    lmrCfg.gridDir = lmrJSONCfg["grid-directory"].str;
    lmrCfg.gridDirectory = lmrCfg.simDir ~ "/" ~ lmrCfg.gridDir;
    lmrCfg.gzipExt = lmrJSONCfg["gzip-extension"].str;
    lmrCfg.blkListFile = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["block-list-filename"].str;
    lmrCfg.vtkDir = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["vtk-output-directory"].str;
    lmrCfg.mpimapFile = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["mpimap-filename"].str;
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
 * Return the flow solution filename for a single block ('blkId') as a string.
 *
 * Authors: RJG
 * Date: 2023-06-27
 */
string flowFilename(int snapshot, int blkId)
{
    string fname = lmrCfg.snapshotDir ~
	"/" ~
	format(lmrCfg.snapshotIdxFmt, snapshot) ~
	"/" ~
	lmrCfg.flowPrefix ~ "-" ~ format(lmrCfg.blkIdxFmt, blkId);
    if (GlobalConfig.flow_format == "gziptext")
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
    if (GlobalConfig.flow_format == "gziptext")
	fname ~= lmrCfg.gzipExt;
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
    if (GlobalConfig.grid_format == "gziptext")
	gname ~= lmrCfg.gzipExt;
    return gname;
}

/** lmrconfig.d
 * Module for configuration of Eilmer program itself.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2022-08-06
 */

module lmrconfig;

import std.format : format;
import std.process : environment;
import std.json;
import json_helper : readJSONfile;

struct LmrCfg {
    immutable string cfgDir;
    immutable string cfgFile;
    immutable string ctrlFile;
    immutable string nkCfgFile;
    immutable string gridDir;
    immutable string blkFmtStr;
    immutable string snapshotDir;
    immutable string snapshotIdxFmtStr;
    immutable string flowDir;
    immutable string zipExt;
    immutable string gzipExt;
    immutable string rawBinExt;
    immutable string gridJobName;
    immutable string flowJobName;
    immutable string blkListFile;
    immutable string vtkDir;
};

LmrCfg lmrCfg;
JSONValue lmrJSONCfg;

static this()
{
    // Read base config into lmrJSONCfg
    readLmrConfig();
    
    // Populate struct with derived config information
    lmrCfg.cfgDir = lmrJSONCfg["config-directory"].str;
    lmrCfg.cfgFile = lmrCfg.cfgDir ~ "/" ~ lmrJSONCfg["config-filename"].str;
    lmrCfg.ctrlFile = lmrCfg.cfgDir ~ "/" ~ lmrJSONCfg["control-filename"].str;
    lmrCfg.nkCfgFile = lmrCfg.cfgDir ~ "/" ~ lmrJSONCfg["newton-krylov-config-filename"].str;
    lmrCfg.gridDir = lmrJSONCfg["grid-directory"].str;
    lmrCfg.blkFmtStr = lmrJSONCfg["block-filename-format"].str;
    lmrCfg.snapshotDir = lmrJSONCfg["snapshot-directory"].str;
    lmrCfg.snapshotIdxFmtStr = lmrJSONCfg["snapshot-index-format"].str;
    lmrCfg.flowDir = lmrJSONCfg["flow-directory"].str;
    lmrCfg.zipExt = lmrJSONCfg["zip-extension"].str;
    lmrCfg.gzipExt = lmrJSONCfg["gziptext-extension"].str;
    lmrCfg.rawBinExt = lmrJSONCfg["rawbinary-extension"].str;
    lmrCfg.gridJobName = lmrJSONCfg["user-supplied-grid-job-name"].str;
    lmrCfg.flowJobName = lmrJSONCfg["user-supplied-flow-job-name"].str;
    lmrCfg.blkListFile = lmrCfg.cfgDir ~ "/" ~ lmrJSONCfg["block-list-filename"].str;
    lmrCfg.vtkDir = lmrJSONCfg["vtk-output-directory"].str;

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
 * Return the grid filename for a single block ('id') as a string.
 *
 * Authors: RJG
 * Date: 2022-08-06
 */
string gridFilenameWithoutExt(int id)
{
    return lmrCfg.gridDir ~ "/" ~ format(lmrCfg.blkFmtStr, id);
}

/**
 * Return the directory name for a snapshot of a flow field.
 *
 * Authors: RJG
 * Date: 2022-08-07
 */
string steadyFlowDirectory(int snapshot)
{
    return lmrCfg.snapshotDir ~ "/" ~ format(lmrCfg.snapshotIdxFmtStr, snapshot) ~ "/" ~ lmrCfg.flowDir;
}

/**
 * Return the flow solution filename for a single block ('id') as a string.
 *
 * Authors: RJG
 * Date: 2022-08-06
 */
string steadyFlowFilename(int snapshot, int blkId)
{
    return steadyFlowDirectory(snapshot) ~ "/" ~ format(lmrCfg.blkFmtStr, blkId) ~ "." ~ lmrCfg.zipExt;
}


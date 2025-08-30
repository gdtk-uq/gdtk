module lmr.lmrbuild;

import std.compiler;
import std.format;
import std.json;

import util.json_helper; 

// Must provide file location with -J option
enum rawBuildInfo = import("buildinfo.json");

// Use enum to allow for compile-time inlining
struct BuildCfg { 

    // We can find these from special variables
    immutable string buildDate;
    enum compilerName = std.compiler.name
            ~ format(" v%d.%d", std.compiler.version_major, std.compiler.version_minor);

    immutable string revisionDate;
    immutable string revisionId;
    immutable string fullRevisionId;

    // Version info
    // We could probably use `debug` and `D_Optimized` versions here.
    // For some reason a `profile` one doesn't exist though :'(
    version (flavour_debug) {
        enum buildFlavour = "debug";
    } else version(flavour_profile) {
        enum buildFlavour = "profile";
    } else version(flavour_fast) {
        enum buildFlavour = "fast";
    } else {
        static assert(0);
    }

    // Are we running in parallel?
    version (mpi_parallel) {
        // Get the parallel implementation from buildinfo.json
        immutable string parallelType;
    } else {
        immutable string parallelType = "shared";
    }

    // Are we using complex numbers?
    version (complex_numbers) {
        enum numberType = "complex";
    } else {
        enum numberType = "real";
    }
}

BuildCfg buildCfg;

static this() {
    JSONValue buildInfo = parseJSON(rawBuildInfo);
    buildCfg.buildDate = getJSONstring(buildInfo, "buildDate", __TIMESTAMP__);
    buildCfg.revisionDate = getJSONstring(buildInfo, "revisionDate", "ERR: No revision date");
    buildCfg.revisionId = getJSONstring(buildInfo, "revisionId", "ERR: No revision id");
    buildCfg.fullRevisionId = getJSONstring(buildInfo, "fullRevisionId", "ERR: No revision id");

    version (mpi_parallel) {
        buildCfg.parallelType = getJSONstring(buildInfo, "parallelType", "ERR: Parallel version incorrectly set");
    }
}

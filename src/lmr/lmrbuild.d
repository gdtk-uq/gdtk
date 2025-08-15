module lmr.lmrbuild;

import std.compiler;

// Use enum to allow for compile-time inlining
struct BuildCfg { 

    // We can find these from special variables
    enum buildDate      = __TIMESTAMP__;
    enum compilerName   = std.compiler.name
            ~ format(" v%d.%d", std.compiler.version_major, std.compiler.version_minor);

    // Git info (commands replaced by shell-script)
    version (none) {
        enum revisionDate   = "{{git log -1 --format=%cd}}";
        enum revisionId     = "{{git rev-parse --short HEAD}}";
        enum fullRevisionId = "{{git rev-parse HEAD}}";
    } else {
        
    }

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
        enum parallelType = "parallel";
    } else {
        enum parallelType = "shared";
    }

    version (complex_numbers) {
        enum numberType = "complex";
    } else {
        enum numberType = "real";
    }
}

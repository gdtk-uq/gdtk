module lmr.lmrbuild;

// Use enum to allow for compile-time inlining
enum buildDate      = "PUT_BUILD_DATE_HERE";
enum buildFlavour   = "PUT_BUILD_FLAVOUR_HERE";
enum compilerName   = "PUT_COMPILER_NAME_HERE";
enum fullRevisionId = "PUT_FULL_REVISION_STRING_HERE";
enum parallelType   = "PUT_PARALLEL_TYPE_HERE";
enum numberType     = "PUT_NUMBER_TYPE_HERE";
enum revisionDate   = "PUT_REVISION_DATE_HERE";
enum revisionId     = "{{git rev-parse --short HEAD}}";

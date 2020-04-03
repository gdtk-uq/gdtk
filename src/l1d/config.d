// config.d
// A place to keep the configuration details.

module config;

final class L1dConfig {
    shared static int verbosity_level = 1;
    // Messages have a hierarchy:
    // 0 : only error messages will be omitted
    // 1 : emit messages that are useful for a long-running job (default)
    // 2 : plus verbose init messages
    // 3 : plus verbose boundary condition messages
    //
    shared static string base_file_name = "job"; // Change this to suit at run time.
}

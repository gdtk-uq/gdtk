module gas.test_gas_calc;

import std.algorithm;
import std.array;
import std.conv;
import std.process;
import std.stdio;

unittest {
    import std.file;

    chdir("./sample-data");
    scope (exit)
        chdir("..");
    
    string goodFileName = "wrapped-gas-model-test.lua";
    int returnCode = wait(spawnProcess(["../gas-calc", goodFileName]));
    assert(returnCode == 0, "Process failed.");
}

// Ensure our tests are isolated
unittest {
    string badFileName = "does-not-exist.lua";
    int returnCode = wait(spawnProcess(["./gas-calc", badFileName]));
    assert(returnCode == 1, "Process succeeded unexpectedly.");
}

// test_runner.d
import core.runtime;
import core.exception;
import std.stdio;
import std.string;
import std.algorithm;
import std.range;
import std.conv;
import std.format;

// Helper functions for adaptive padding
string centerPadded(string text, char padChar, int totalWidth = 79) {
    int textLen = cast(int)text.length;
    if (textLen >= totalWidth) return text;
    
    int padding = totalWidth - textLen;
    int leftPad = padding / 2;
    int rightPad = padding - leftPad;
    
    return padChar.repeat(leftPad).to!string ~ text ~ padChar.repeat(rightPad).to!string;
}

string leftPadded(string text, char padChar, int totalWidth = 79) {
    int textLen = cast(int)text.length;
    if (textLen >= totalWidth) return text;
    
    int padding = totalWidth - textLen;
    return padChar.repeat(padding).to!string ~ text;
}

// ANSI color codes
immutable string RED = "\033[31m";
immutable string GREEN = "\033[32m";
immutable string YELLOW = "\033[33m";
immutable string BLUE = "\033[34m";
immutable string MAGENTA = "\033[35m";
immutable string CYAN = "\033[36m";
immutable string WHITE = "\033[37m";
immutable string BOLD = "\033[1m";
immutable string RESET = "\033[0m";

struct TestStats {
    size_t passed;
    size_t failed;
    string[] failures;
    double duration;
}

UnitTestResult customUnitTester()
{
    import std.datetime.stopwatch;
    import core.runtime : UnitTestResult;
    
    auto sw = StopWatch(AutoStart.yes);
    TestStats stats;
    
    // Header like pytest
    writefln(BOLD ~ centerPadded(" test session starts ", '=') ~ RESET);
    
    // Collect modules with unit tests
    ModuleInfo*[] modules;
    foreach (m; ModuleInfo) {
        if (m !is null && m.unitTest !is null) {
            modules ~= m;
        }
    }
    
    if (modules.length == 0) {
        writeln("No tests found.");
        return UnitTestResult(0, 0, false, true);
    }
    
    writefln("collected %d test module%s", modules.length, modules.length == 1 ? "" : "s");
    writeln();
    
    // Progress bar style output
    write("Running tests: ");
    stdout.flush();
    
    foreach (i, m; modules) {
        string moduleName = m.name;
        
        // Shorten module names for cleaner display
        if (moduleName.canFind('.')) {
            auto parts = moduleName.split('.');
            moduleName = parts[$ - 1]; // Just the last part
        }
        
        try {
            m.unitTest()();  // Call the function pointer correctly
            write(GREEN ~ "." ~ RESET);
            stats.passed++;
        } catch (AssertError e) {
            write(RED ~ "F" ~ RESET);
            stats.failed++;
            
            // Store failure info for later display
            auto failureMsg = format("%s::%s\n    %s:%s - %s", 
                                   moduleName, 
                                   "unittest", 
                                   e.file.split('/').back.split('\\').back, // Just filename
                                   e.line, 
                                   e.msg.length ? e.msg : "assertion failed");
            stats.failures ~= failureMsg;
        } catch (Throwable t) {
            write(YELLOW ~ "E" ~ RESET);
            stats.failed++;
            
            auto failureMsg = format("%s::%s\n    ERROR: %s", 
                                   moduleName, 
                                   "unittest", 
                                   t.msg);
            stats.failures ~= failureMsg;
        }
        
        stdout.flush();
    }
    
    sw.stop();
    stats.duration = sw.peek().total!"msecs" / 1000.0;
    
    writeln();
    writeln();
    
    // Show failures if any
    if (stats.failed > 0) {
        writefln(RED ~ BOLD ~ centerPadded(" FAILURES ", '=') ~ RESET);
        
        foreach (i, failure; stats.failures) {
            if (i > 0) writeln();
            writefln(RED ~ centerPadded(" " ~ failure.split('\n')[0] ~ " ", '_') ~ RESET);
            writeln(failure.split('\n')[1..$].join("\n"));
        }
        writeln();
    }
    
    // Summary like pytest
    string summaryColor = stats.failed > 0 ? RED : GREEN;
    string status = stats.failed > 0 ? "FAILED" : "PASSED";
    
    write(summaryColor ~ BOLD);
    
    if (stats.failed > 0) {
        string summaryText = format(" %s, %d failed, %d passed in %.2fs ", 
                                   status, stats.failed, stats.passed, stats.duration);
        writef(summaryColor ~ BOLD ~ centerPadded(summaryText, '=') ~ RESET);
    } else {
        string summaryText = format(" %d passed in %.2fs ", stats.passed, stats.duration);
        writef(summaryColor ~ BOLD ~ centerPadded(summaryText, '=') ~ RESET);
    }
    
    writeln(RESET);
    
    return UnitTestResult(
        stats.passed + stats.failed,  // executed
        stats.passed,                 // passed
        false,                       // summarize - suppress default message
        stats.failed == 0            // result
    );
}

static this()
{
    Runtime.extendedModuleUnitTester = &customUnitTester;
}

void main()
{
    writeln("Test discovery and execution handled by custom runner...");
}

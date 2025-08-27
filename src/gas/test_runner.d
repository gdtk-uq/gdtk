// test_runner.d
import core.exception;
import core.runtime;

import std.algorithm;
import std.conv;
import std.format;
import std.getopt;
import std.range;
import std.stdio;
import std.string;

import gas;

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

    string[] args = Runtime.args;
    string rootModule = "";
    string[] excludePatterns;
    auto opts = getopt(
        args,
        "module", &rootModule,
        "exclude", &excludePatterns,
    );

    auto sw = StopWatch(AutoStart.yes);
    TestStats stats;
    
    // Header like pytest
    writefln(BOLD ~ centerPadded(" test session starts ", '=') ~ RESET);
    
    // Collect modules with unit tests
    ModuleInfo*[] modules;
    foreach (m; ModuleInfo) {
        if (!(m.name).startsWith(rootModule)) {
            continue;
        }
        if (m !is null && m.unitTest !is null) {
            modules ~= m;
        }
    }

    // Sort the modules so they are easier to look through
    modules.sort!((a, b) => a.name < b.name);
    
    if (modules.length == 0) {
        writeln("No tests found.");
        return UnitTestResult(0, 0, false, true);
    }
    
    writefln("collected %d test module%s", modules.length, modules.length == 1 ? "" : "s");
    writeln();
    
    // Verbose pytest-style output
    foreach (i, m; modules) {
        string testName = m.name ~ "::unittest";
        write(testName ~ " ... ");
        stdout.flush();
        
        try {
            auto testFunc = m.unitTest();
            testFunc();
            writeln(GREEN ~ "PASSED" ~ RESET);
            stats.passed++;
        } catch (AssertError e) {
            writeln(RED ~ "FAILED" ~ RESET);
            stats.failed++;
            
            // Store failure info for later display
            auto failureMsg = format("%s\n    %s:L%s - %s", 
                                   testName,
                                   e.file.split('/').back.split('\\').back, // Just filename
                                   e.line, 
                                   e.msg.length ? e.msg : "assertion failed");
            stats.failures ~= failureMsg;
        } catch (Throwable t) {
            writeln(YELLOW ~ "ERROR" ~ RESET);
            stats.failed++;
            
            auto failureMsg = format("%s\n    ERROR: %s", 
                                   testName,
                                   t.msg);
            stats.failures ~= failureMsg;
        }
    }
    
    sw.stop();
    stats.duration = sw.peek().total!"msecs" / 1000.0;
    
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
    
    if (stats.failed > 0) {
        string summaryText = format(" %s, %d failed, %d passed in %.2fs ", 
                                   status, stats.failed, stats.passed, stats.duration);
        writeln(summaryColor ~ BOLD ~ centerPadded(summaryText, '=') ~ RESET);
    } else {
        string summaryText = format(" %d passed in %.2fs ", stats.passed, stats.duration);
        writeln(summaryColor ~ BOLD ~ centerPadded(summaryText, '=') ~ RESET);
    }
    
    return UnitTestResult(
        executed: stats.passed + stats.failed,  // executed
        passed: stats.passed,                 // passed
        summarize: false,                       // summarize - suppress default message
        runMain: true,
        // stats.failed == 0            // result
    );
}

static this()
{
    Runtime.extendedModuleUnitTester = &customUnitTester;
}

void main() {}

// test_runner.d
module util.test_runner;

import core.exception : AssertError;
import core.runtime : Runtime, UnitTestResult;

import std.ascii;
import std.algorithm;
import std.conv;
import std.format;
import std.getopt;
import std.path : baseName;
import std.range;
import std.regex;
import std.stdio;
import std.string;

immutable OUTPUT_WIDTH = 79;

// ANSI color codes
immutable string RED = "\x1b[31m";
immutable string GREEN = "\x1b[32m";
immutable string YELLOW = "\x1b[33m";
immutable string BLUE = "\x1b[34m";
immutable string MAGENTA = "\x1b[35m";
immutable string CYAN = "\x1b[36m";
immutable string WHITE = "\x1b[37m";
immutable string BOLD = "\x1b[1m";
immutable string RESET = "\x1b[0m";

// Allows us to detect skips, and differentiate them from assert errors
// However, with a standard test runner, it will still throw errors, meaning
// they can be caught.
class SkipError : AssertError {
    @safe pure nothrow @nogc this(string file, size_t line) {
        this(cast(Throwable) null, file, line);
    }

    @safe pure nothrow @nogc this(Throwable next, string file = __FILE__, size_t line = __LINE__) {
        this("Skipped test", file, line, next);
    }

    @safe pure nothrow @nogc this(string msg, string file = __FILE__, size_t line = __LINE__, Throwable next = null) {
        super(msg, file, line, next);
    }
}

void skip(string message = "Skipped test") {
    throw new SkipError(message);
}

size_t visualWidth(string text) {
    static auto ansiRegex = ctRegex!(r"\x1b\[[0-9;]*[A-Za-z]");
    return replaceAll(text, ansiRegex, "").length;
}

/**
 * Drop-in replacement for std.string.center that accounts for ANSI escape sequences
 */
string center(string text, size_t width, dchar fillChar = ' ') {
    size_t ansiOverhead = text.length - visualWidth(text);
    return std.string.center(text, width + ansiOverhead, fillChar);
}

/**
 * Drop-in replacement for std.string.leftJustify that accounts for ANSI escape sequences
 */
string leftJustify(string text, size_t width, dchar fillChar = ' ') {
    size_t ansiOverhead = text.length - visualWidth(text);
    return std.string.leftJustify(text, width + ansiOverhead, fillChar);
}

enum TestResult {
    PASSED = 0,
    FAILED,
    SKIPPED,
    ERROR,
}

string getColour(TestResult r) {
    final switch (r) {
    case TestResult.PASSED:
        return GREEN;
    case TestResult.FAILED:
        return RED;
    case TestResult.SKIPPED, TestResult.ERROR:
        return YELLOW;
    }
}

struct TestInfo {
    string name;
    TestResult result;
    string context = "";
    string message = "";

    string toString() const {
        if (result == TestResult.PASSED) {
            return "";
        }

        string output = result.getColour();
        output ~= center(" " ~ name ~ " ", OUTPUT_WIDTH, '_') ~ RESET ~ "\n";
        output ~= context ~ (context.length ? "\n" : "");
        output ~= "  " ~ BOLD ~ message ~ RESET ~ "\n";

        return output;
    }
}

struct Summary {
    size_t passed;
    size_t failed;
    size_t skipped;
    double duration;

    string toString() const {
        string summaryText = "";
        string summaryColour = failed > 0 ? RED : GREEN;

        if (failed > 0) {
            summaryText ~= BOLD ~ format("%d failed", failed) ~ RESET ~ ", ";
        }
        summaryText ~= GREEN ~ format("%d passed", passed) ~ summaryColour;
        if (skipped > 0) {
            summaryText ~= RESET ~ ", " ~ YELLOW ~ format("%d skipped", skipped) ~ summaryColour;
        }
        summaryText ~= format(" in %.2fs", duration);
        return summaryColour ~ center(" " ~ summaryText ~ " ", OUTPUT_WIDTH, '=') ~ RESET;
    }
}

TestInfo captureTest(void function() testFunc) {
    // Set up the redirect from stdout to File
    auto oldStdout = stdout;
    scope (exit) // Ensure we put everything back when we're done
        stdout = oldStdout;

    stdout = File.tmpfile();

    TestInfo testInfo;

    try {
        testFunc();
        testInfo.result = TestResult.PASSED;
    } catch (SkipError e) {
        testInfo.message = format("%s:L%s - %s", baseName(e.file), e.line, e.msg.length ? e.msg : "test skipped");
        testInfo.result = TestResult.SKIPPED;
    } catch (AssertError e) {
        // Store failure info for later display
        testInfo.message = format("%s:L%s - %s", baseName(e.file), e.line,
                e.msg.length ? e.msg : "assertion failed");
        testInfo.result = TestResult.FAILED;
    } catch (Throwable t) {
        testInfo.message ~= format("%s:L%s - ERROR: %s", baseName(t.file), t.line, t.msg);
        testInfo.result = TestResult.ERROR;
    } finally {
        stdout.flush();
        stdout.rewind();
        testInfo.context = stdout.byLine().map!text.join("\n");
    }

    return testInfo;
}

UnitTestResult customUnitTester() {
    import std.datetime.stopwatch;

    string[] args = Runtime.args;
    string rootModule = "";
    string[] excludePatterns;
    auto opts = getopt(args, "module", &rootModule, "exclude", &excludePatterns);

    // Header like pytest
    writefln(BOLD ~ center(" test session starts ", OUTPUT_WIDTH, '=') ~ RESET);

    ModuleInfo*[] allModules;
    foreach (m; ModuleInfo) {
        allModules ~= m;
    }

    // Sort and filter the modules
    auto modules = allModules.sort!((a, b) => a.name < b.name)
        .filter!(m => m.unitTest !is null)
        .filter!(m => (m.name).startsWith(rootModule))
        .array;

    if (modules.length == 0) {
        writeln("No tests found.");
        return UnitTestResult(0, 0, false, true);
    }

    writefln("collected %d test module%s", modules.length, modules.length == 1 ? "" : "s");
    writeln();

    auto sw = StopWatch(AutoStart.yes);
    TestInfo[] allTests;

    // Verbose pytest-style output
    foreach (i, m; modules) {
        string testName = m.name ~ "::unittest";
        write(testName ~ " ... ");
        stdout.flush();

        auto testFunc = m.unitTest();
        auto testInfo = captureTest(testFunc);
        testInfo.name = testName;

        write("\b\b\b\b\b "); // Erase the string " ... "
        write(testInfo.result.getColour());
        write(testInfo.result);
        writeln(RESET);

        allTests ~= testInfo;
    }

    sw.stop();
    auto failedTests = allTests.filter!(t => (t.result == TestResult.ERROR) || (t.result == TestResult.FAILED)).array();
    auto numSkippedTests = allTests.filter!(t => t.result == TestResult.SKIPPED).count();
    auto numPassedTests = allTests.filter!(t => t.result == TestResult.PASSED).count();

    writeln();

    // Show failures if any
    if (failedTests.length > 0) {
        writefln(BOLD ~ center(" FAILURES ", OUTPUT_WIDTH, '=') ~ RESET);

        foreach (failure; failedTests) {
            writeln(failure);
        }
    }

    Summary summary = {
        duration: sw.peek().total!"msecs" / 1000.0,
        passed: numPassedTests,
        failed: failedTests.length,
        skipped: numSkippedTests
    };
    writeln(summary);

    UnitTestResult result = {
        executed: summary.passed + summary.failed, passed: summary.passed,
        summarize: false, runMain: true
    };

    return result;
}

version (unittest) {
    static this() {
        Runtime.extendedModuleUnitTester = &customUnitTester;
    }

    void main() {
    }
}

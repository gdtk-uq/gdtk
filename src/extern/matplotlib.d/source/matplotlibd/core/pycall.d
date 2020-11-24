module matplotlibd.core.pycall;


void call(string py_script) {
    import std.process: environment, pipeProcess, wait, Redirect;

    auto py_path =  environment.get("MATPLOTLIB_D_PYTHON", "python3");
    auto pipes = pipeProcess(py_path, Redirect.stdin | Redirect.stderr);

    pipes.stdin.writeln(py_script);
    pipes.stdin.writeln("exit()");
    pipes.stdin.close();

    auto result = wait(pipes.pid);

    if (result != 0) {
        string error;
        foreach (line; pipes.stderr.byLine)
            error ~= line ~ "\n";
        throw new Exception("\n\nERROR occurred in Python:\n" ~ error);
    }
}

unittest {
    call("print(\"Passing!\")\n");
}

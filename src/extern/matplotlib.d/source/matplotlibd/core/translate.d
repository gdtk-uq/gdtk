module matplotlibd.core.translate;

alias immutable bool PyBool;
alias immutable (void*) PyNone;

PyBool False = false;
PyBool True = true;
PyNone None = null;


string d2py(T)(T v) {
    import std.format: format;
    static if (is(typeof(v) : PyNone))
        return "None";

    else static if (is(typeof(v) : bool))
        return v ? "True" : "False";        

    else static if (is(typeof(v) : string))
        return format("\"%s\"", v);

    else
        return format("%s", v);
}

unittest {
    import std.range: iota;
    assert(d2py(None) == "None");
    assert(d2py(null) == "None");
    assert(d2py(True) == "True");
    assert(d2py(true) == "True");
    assert(d2py(False) == "False");
    assert(d2py(false) == "False");
    assert(d2py("Hello!") == "\"Hello!\"");
    assert(d2py(5.iota) == "[0, 1, 2, 3, 4]");
}


string parseArgs(Args)(Args args) {
    static if (is(typeof(args.keys) : string[])) {
        string parsed;
        foreach(key; args.byKey)
            parsed ~= key ~ "=" ~  d2py(args[key]) ~ ",";
    }
    else
        string parsed =  d2py(args) ~ ",";
    return parsed;
}

unittest {
    import std.range: iota;
    assert(parseArgs(5) == "5,");
    assert(parseArgs(5.iota) == "[0, 1, 2, 3, 4],");
    assert(parseArgs(["test": 5]) == "test=5,");
    assert(parseArgs(["test": "test"]) == "test=\"test\",");
    assert(parseArgs(["test": 5.iota]) == "test=[0, 1, 2, 3, 4],");
    assert(parseArgs(["test": false]) == "test=False,");
    assert(parseArgs(["test": False]) == "test=False,");
}

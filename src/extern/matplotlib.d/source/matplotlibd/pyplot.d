module matplotlibd.pyplot;
import matplotlibd.core.pycall;
import matplotlibd.core.translate;

private:

string py_script = "import matplotlib.pyplot as plt\n";

immutable string plt_funcs = (){
    import std.string: splitLines;

    string[] method_names = import("pyplot_functions.txt").splitLines;
    string plt_funcs;

    foreach(name; method_names) {
        plt_funcs ~=
            "void " ~ name ~ "(T...)(T a)" ~
            "{import std.format: format;" ~
            "string p;if(a.length>0){foreach(i;a){p~=parseArgs(i);}" ~
            "p = p[0..$-1];}py_script~=format(\"plt."~ name ~ "(%s)\n\",p);";

        if (name == "show" || name == "savefig")
            plt_funcs ~= "call(py_script);}\n";
        else
            plt_funcs ~= "}\n";
    }

    return plt_funcs;
}();


public:

import matplotlibd.core.translate: False, True, None;


void clear() {
    py_script = "import matplotlib.pyplot as plt\n";
}

mixin(plt_funcs);

unittest {
    import std.string;
    auto script = py_script ~ "plt.plot([1, 2],[2, 4],\"r-\",lw=2)\n";
    plot([1, 2], [2, 4], "r-", ["lw": 2]);
    assert(script == py_script);
    clear();
    assert(py_script == "import matplotlib.pyplot as plt\n");
}

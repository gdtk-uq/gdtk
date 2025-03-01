module force_tuple_instantiation;

import std.typecons;

static this() {
    // Force the instantiation of Tuple!(int, double).opEquals by invoking it.
    Tuple!(int, double) a, b;
    a.opEquals(b);
}


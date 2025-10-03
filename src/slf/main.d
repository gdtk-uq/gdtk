/*
    main.d: a basic executable entry-point to the slf calculator
*/

import std.string;

import slf.core : Flame;

int main(string[] args)
{
    int exitFlag = 0;
    string name = "slf";
    if (args.length>1) name = args[1];
    name = name.chomp(".yaml");

    auto flame = new Flame(name);
    flame.set_initial_condition();

    exitFlag = flame.run();

    flame.save_solution();
    flame.save_log();

    return exitFlag;
} // end main()

module lmr.commands.command;

import core.stdc.stdlib : system;
import std.string;

enum LmrCmdType {
    user,
    dev
};

struct Command
{
    /// pointer to function for command action
    int function(string[]) main;

    string description;
    string shortDescription;
    string helpMsg;
    LmrCmdType type = LmrCmdType.user; // by default


}

class UserInputError : Error {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

string shellCommand(string[] args)
{
    string str = args[0] ~ "-" ~ args[1];
    foreach (s; args[2 .. $]) {
        str ~= " " ~ s;
    }
    return str;

}

int callShellCommand(string[] args)
{
    auto shellCmd = shellCommand(args);
    return system(shellCmd.toStringz);
}

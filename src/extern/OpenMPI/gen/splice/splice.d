import std.file : readText;
import std.stdio, std.algorithm, std.range, std.string;

void main(string[] args)
{
    auto pack = readText(args[1]);
    auto cfg = readText(args[2]).splitLines;

    auto r0 = pack.findSplitAfter("// BEGIN AUTO\n");
    auto r1 = r0[1].findSplitAfter("// BEGIN AUTO\n");
    auto pre = r0[0];
    auto mid = r1[0];
    auto post = r1[1];
    chain(pre,
            cfg.filter!(s => !s.canFind("OMPI_OFFSET_DATATYPE")).joiner("\n"),
            mid,
            cfg.filter!(s => s.canFind("OMPI_OFFSET_DATATYPE")).joiner("\n"),
            post)
        .write();
}

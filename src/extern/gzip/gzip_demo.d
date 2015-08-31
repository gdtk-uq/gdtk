import gzip;
import std.stdio;

void main() {
    auto byLine = new GzipByLine("test.txt.gz");
    foreach(line; byLine)
	writeln(line);

    auto gzipOutFile = new GzipOut("testout.gz");
    gzipOutFile.compress("bla bla bla");
    gzipOutFile.finish();
}


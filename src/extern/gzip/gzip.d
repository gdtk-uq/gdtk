// gzip.d
// Python-inspired handling of gzip files.
// Author: Stephan Schiffels (on forum.dlang.org)
// Date: 2014-feb-20
//
// cfcfd3 version: 2014-07-21, PJ

/+
Hi Kamil,
I am glad someone has the exact same problem as I had. I actually 
solved this, inspired by the python API you quoted above. I wrote 
these classes:
GzipInputRange, GzipByLine, and GzipOut.
Here is how I can now use them:

_____________________

import gzip;
import std.stdio;

void main() {
auto byLine = new GzipByLine("test.gz");
foreach(line; byLine)
   writeln(line);

auto gzipOutFile = new GzipOut("testout.gz");
gzipOutFile.compress("bla bla bla");
gzipOutFile.finish();
}

That is all quite convenient and I was wondering whether 
something like that would be useful even in Phobos. But it's 
clear that for phobos things would involve a lot more work to 
comply with the requirements. This so far simply served my needs 
and is not as generic as it could be:

Here is the code:

+/

module gzip;

import std.zlib;
import std.stdio;
import std.range;
import std.traits;

class GzipInputRange {
  UnCompress uncompressObj;
  File f;
  auto CHUNKSIZE = 0x4000;
  ReturnType!(f.byChunk) chunkRange;
  bool exhausted;
  char[] uncompressedBuffer;
  size_t bufferIndex;

  this(string filename) {
    f = File(filename, "r");
    chunkRange = f.byChunk(CHUNKSIZE);
    uncompressObj = new UnCompress();
    load();
  }

  void load() {
    if(!chunkRange.empty) {
      auto raw = chunkRange.front.dup;
      chunkRange.popFront();
      uncompressedBuffer = cast(char[])uncompressObj.uncompress(raw);
      bufferIndex = 0;
    }
    else {
      if(!exhausted) {
        uncompressedBuffer = cast(char[])uncompressObj.flush();
        exhausted = true;
        bufferIndex = 0;
      }
      else
        uncompressedBuffer.length = 0;
    }
  }

  @property char front() {
    return uncompressedBuffer[bufferIndex];
  }

  void popFront() {
    bufferIndex += 1;
    if(bufferIndex >= uncompressedBuffer.length) {
      load();
      bufferIndex = 0;
    }
  }

  @property bool empty() {
    return uncompressedBuffer.length == 0;
  }
}

class GzipByLine {
  GzipInputRange range;
  char[] buf;

  this(string filename) {
    this.range = new GzipInputRange(filename);
    popFront();
  }

  @property bool empty() {
    return buf.length == 0;
  }

  void popFront() {
    buf.length = 0;
    while(!range.empty && range.front != '\n') {
      buf ~= range.front;
      range.popFront();
    }
    range.popFront();
  }

  string front() {
    return buf.idup;
  }
}

class GzipOut {
  Compress compressObj;
  File f;

  this(string filename) {
    f = File(filename, "w");
    compressObj = new Compress(HeaderFormat.gzip);
  }

  void compress(string s) {
    auto compressed = compressObj.compress(s.dup);
    f.rawWrite(compressed);
  }

  void finish() {
    auto compressed = compressObj.flush();
    f.rawWrite(compressed);
  }
}

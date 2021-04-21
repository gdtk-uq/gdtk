/*
  D bindings for gperftools(Google Performance Tools).
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module gperftools_d.heap_profiler;

extern (C):

void HeapProfilerStart (const(char)* prefix);

int IsHeapProfilerRunning ();

void HeapProfilerStop ();

void HeapProfilerDump (const(char)* reason);

char* GetHeapProfile ();

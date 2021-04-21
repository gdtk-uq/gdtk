/*
  D bindings for gperftools(Google Performance Tools).
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module gperftools_d.stacktrace;

extern (C):

int GetStackFrames (void** result, int* sizes, int max_depth, int skip_count);

int GetStackFramesWithContext (
    void** result,
    int* sizes,
    int max_depth,
    int skip_count,
    const(void)* uc);

int GetStackTrace (void** result, int max_depth, int skip_count);

int GetStackTraceWithContext (
    void** result,
    int max_depth,
    int skip_count,
    const(void)* uc);

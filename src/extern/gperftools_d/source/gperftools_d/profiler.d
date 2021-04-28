/*
  D bindings for gperftools(Google Performance Tools).
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module gperftools_d.profiler;

import core.stdc.time;

extern (C):


struct ProfilerOptions
{
    int function (void* arg) filter_in_thread;
    void* filter_in_thread_arg;
}

int ProfilerStart ();


int ProfilerStartWithOptions (
    const(char)* fname,
    const(ProfilerOptions)* options);

void ProfilerStop ();

void ProfilerFlush ();

void ProfilerEnable ();
void ProfilerDisable ();

int ProfilingIsEnabledForAllThreads ();

void ProfilerRegisterThread ();

struct ProfilerState
{
    int enabled;
    time_t start_time;
    char[1024] profile_name;
    int samples_gathered;
}

void ProfilerGetCurrentState (ProfilerState* state);


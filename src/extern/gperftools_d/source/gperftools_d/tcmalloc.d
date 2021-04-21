/*
  D bindings for gperftools(Google Performance Tools).
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module gperftools_d.tc_malloc;

extern (C):

enum TC_VERSION_MAJOR = 2;
enum TC_VERSION_MINOR = 4;
enum TC_VERSION_PATCH = "";
enum TC_VERSION_STRING = "gperftools 2.4";

const(char)* tc_version (int* major, int* minor, const(char*)* patch);

void* tc_malloc (size_t size);
void* tc_malloc_skip_new_handler (size_t size);
void tc_free (void* ptr);
void* tc_realloc (void* ptr, size_t size);
void* tc_calloc (size_t nmemb, size_t size);
void tc_cfree (void* ptr);

void* tc_memalign (size_t __alignment, size_t __size);
int tc_posix_memalign (void** ptr, size_t align_, size_t size);
void* tc_valloc (size_t __size);
void* tc_pvalloc (size_t __size);

void tc_malloc_stats ();
int tc_mallopt (int cmd, int value);

size_t tc_malloc_size (void* ptr);

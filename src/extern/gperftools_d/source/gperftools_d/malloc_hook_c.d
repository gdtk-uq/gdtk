/*
  D bindings for gperftools(Google Performance Tools).
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module gperftools_d.malloc_hook_c;

extern (C):

int MallocHook_GetCallerStackTrace (
    void** result,
    int max_depth,
    int skip_count);


alias MallocHook_NewHook = void function (const(void)* ptr, size_t size);
int MallocHook_AddNewHook (MallocHook_NewHook hook);
int MallocHook_RemoveNewHook (MallocHook_NewHook hook);

alias MallocHook_DeleteHook = void function (const(void)* ptr);
int MallocHook_AddDeleteHook (MallocHook_DeleteHook hook);
int MallocHook_RemoveDeleteHook (MallocHook_DeleteHook hook);

alias MallocHook_PreMmapHook = void function (const(void)* start, size_t size, int protection, int flags, int fd, int offset);
int MallocHook_AddPreMmapHook (MallocHook_PreMmapHook hook);
int MallocHook_RemovePreMmapHook (MallocHook_PreMmapHook hook);

alias MallocHook_MmapHook = void function (const(void)* result, const(void)* start, size_t size, int protection, int flags, int fd, int offset);
int MallocHook_AddMmapHook (MallocHook_MmapHook hook);
int MallocHook_RemoveMmapHook (MallocHook_MmapHook hook);

alias MallocHook_MmapReplacement = int function (const(void)* start, size_t size, int protection, int flags, int fd, int offset, void** result);
int MallocHook_SetMmapReplacement (MallocHook_MmapReplacement hook);
int MallocHook_RemoveMmapReplacement (MallocHook_MmapReplacement hook);

alias MallocHook_MunmapHook = void function (const(void)* ptr, size_t size);
int MallocHook_AddMunmapHook (MallocHook_MunmapHook hook);
int MallocHook_RemoveMunmapHook (MallocHook_MunmapHook hook);

alias MallocHook_MunmapReplacement = int function (const(void)* ptr, size_t size, int* result);
int MallocHook_SetMunmapReplacement (MallocHook_MunmapReplacement hook);
int MallocHook_RemoveMunmapReplacement (MallocHook_MunmapReplacement hook);

alias MallocHook_MremapHook = void function (const(void)* result, const(void)* old_addr, size_t old_size, size_t new_size, int flags, const(void)* new_addr);
int MallocHook_AddMremapHook (MallocHook_MremapHook hook);
int MallocHook_RemoveMremapHook (MallocHook_MremapHook hook);

alias MallocHook_PreSbrkHook = void function (ptrdiff_t increment);
int MallocHook_AddPreSbrkHook (MallocHook_PreSbrkHook hook);
int MallocHook_RemovePreSbrkHook (MallocHook_PreSbrkHook hook);

alias MallocHook_SbrkHook = void function (const(void)* result, ptrdiff_t increment);
int MallocHook_AddSbrkHook (MallocHook_SbrkHook hook);
int MallocHook_RemoveSbrkHook (MallocHook_SbrkHook hook);

MallocHook_NewHook MallocHook_SetNewHook (MallocHook_NewHook hook);
MallocHook_DeleteHook MallocHook_SetDeleteHook (MallocHook_DeleteHook hook);
MallocHook_PreMmapHook MallocHook_SetPreMmapHook (MallocHook_PreMmapHook hook);
MallocHook_MmapHook MallocHook_SetMmapHook (MallocHook_MmapHook hook);
MallocHook_MunmapHook MallocHook_SetMunmapHook (MallocHook_MunmapHook hook);
MallocHook_MremapHook MallocHook_SetMremapHook (MallocHook_MremapHook hook);
MallocHook_PreSbrkHook MallocHook_SetPreSbrkHook (MallocHook_PreSbrkHook hook);
MallocHook_SbrkHook MallocHook_SetSbrkHook (MallocHook_SbrkHook hook);

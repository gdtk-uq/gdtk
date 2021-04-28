/*
  D bindings for gperftools(Google Performance Tools).
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module gperftools_d.malloc_extension_c;

extern (C):

enum kMallocExtensionHistogramSize = 64;

int MallocExtension_VerifyAllMemory ();
int MallocExtension_VerifyNewMemory (const(void)* p);
int MallocExtension_VerifyArrayNewMemory (const(void)* p);
int MallocExtension_VerifyMallocMemory (const(void)* p);
int MallocExtension_MallocMemoryStats (
    int* blocks,
    size_t* total,
    int[kMallocExtensionHistogramSize] histogram);
void MallocExtension_GetStats (char* buffer, int buffer_length);

int MallocExtension_GetNumericProperty (const(char)* property, size_t* value);
int MallocExtension_SetNumericProperty (const(char)* property, size_t value);
void MallocExtension_MarkThreadIdle ();
void MallocExtension_MarkThreadBusy ();
void MallocExtension_ReleaseToSystem (size_t num_bytes);
void MallocExtension_ReleaseFreeMemory ();
size_t MallocExtension_GetEstimatedAllocatedSize (size_t size);
size_t MallocExtension_GetAllocatedSize (const(void)* p);

enum MallocExtension_Ownership
{
    MallocExtension_kUnknownOwnership = 0,
    MallocExtension_kOwned = 1,
    MallocExtension_kNotOwned = 2
}

MallocExtension_Ownership MallocExtension_GetOwnership (const(void)* p);


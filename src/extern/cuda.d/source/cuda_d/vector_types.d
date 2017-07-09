/*
  D bindings for CUDA.
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module cuda_d.vector_types;

import core.stdc.config;

extern (C):

/*DEVICE_BUILTIN*/
struct char1
{
    byte x;
}

/*DEVICE_BUILTIN*/
struct uchar1
{
    ubyte x;
}

/*DEVICE_BUILTIN*/
struct char2
{
    byte x;
    byte y;
}

/*DEVICE_BUILTIN*/
struct uchar2
{
    ubyte x;
    ubyte y;
}

/*DEVICE_BUILTIN*/
struct char3
{
    byte x;
    byte y;
    byte z;
}

/*DEVICE_BUILTIN*/
struct uchar3
{
    ubyte x;
    ubyte y;
    ubyte z;
}

/*DEVICE_BUILTIN*/
struct char4
{
    byte x;
    byte y;
    byte z;
    byte w;
}

/*DEVICE_BUILTIN*/
struct uchar4
{
    ubyte x;
    ubyte y;
    ubyte z;
    ubyte w;
}

/*DEVICE_BUILTIN*/
struct short1
{
    short x;
}

/*DEVICE_BUILTIN*/
struct ushort1
{
    ushort x;
}

/*DEVICE_BUILTIN*/
struct short2
{
    short x;
    short y;
}

/*DEVICE_BUILTIN*/
struct ushort2
{
    ushort x;
    ushort y;
}

/*DEVICE_BUILTIN*/
struct short3
{
    short x;
    short y;
    short z;
}

/*DEVICE_BUILTIN*/
struct ushort3
{
    ushort x;
    ushort y;
    ushort z;
}

/*DEVICE_BUILTIN*/
struct short4
{
    short x;
    short y;
    short z;
    short w;
}

/*DEVICE_BUILTIN*/
struct ushort4
{
    ushort x;
    ushort y;
    ushort z;
    ushort w;
}

/*DEVICE_BUILTIN*/
struct int1
{
    int x;
}

/*DEVICE_BUILTIN*/
struct uint1
{
    uint x;
}

/*DEVICE_BUILTIN*/
struct int2
{
    int x;
    int y;
}

/*DEVICE_BUILTIN*/
struct uint2
{
    uint x;
    uint y;
}

/*DEVICE_BUILTIN*/
struct int3
{
    int x;
    int y;
    int z;
}

/*DEVICE_BUILTIN*/
struct uint3
{
    uint x;
    uint y;
    uint z;
}

/*DEVICE_BUILTIN*/
struct int4
{
    int x;
    int y;
    int z;
    int w;
}

/*DEVICE_BUILTIN*/
struct uint4
{
    uint x;
    uint y;
    uint z;
    uint w;
}

/*DEVICE_BUILTIN*/
struct long1
{
    c_long x;
}

/*DEVICE_BUILTIN*/
struct ulong1
{
    c_ulong x;
}

/*DEVICE_BUILTIN*/
struct long2
{
    c_long x;
    c_long y;
}

/*DEVICE_BUILTIN*/
struct ulong2
{
    c_ulong x;
    c_ulong y;
}

/*DEVICE_BUILTIN*/
struct long3
{
    c_long x;
    c_long y;
    c_long z;
}

/*DEVICE_BUILTIN*/
struct ulong3
{
    c_ulong x;
    c_ulong y;
    c_ulong z;
}

/*DEVICE_BUILTIN*/
struct long4
{
    c_long x;
    c_long y;
    c_long z;
    c_long w;
}

/*DEVICE_BUILTIN*/
struct ulong4
{
    c_ulong x;
    c_ulong y;
    c_ulong z;
    c_ulong w;
}

/*DEVICE_BUILTIN*/
struct float1
{
    float x;
}

/*DEVICE_BUILTIN*/
struct float2_
{
    float x;
    float y;
}

/*DEVICE_BUILTIN*/
struct float3
{
    float x;
    float y;
    float z;
}

/*DEVICE_BUILTIN*/
struct float4
{
    float x;
    float y;
    float z;
    float w;
}

/*DEVICE_BUILTIN*/
struct double2_
{
    double x;
    double y;
}

/*DEVICE_BUILTIN*/

/*DEVICE_BUILTIN*/
struct dim3
{
    uint x;
    uint y;
    uint z;
}
// This is the REAL "hello world" for CUDA!
// It takes the string "Hello ", prints it, then passes it to CUDA with an array
// of offsets. Then the offsets are added in parallel to produce the string "World!"
// By Ingemar Ragnemalm 2010
// modified by Kyle Damm 2017 for use with D

extern (C) void hello(char *ad, int *bd, int blocksize);

import cuda_d.cuda;
import cuda_d.cuda_runtime_api;
import cuda_d.cublas_api;
import std.stdio;

const int N = 7;
const int blocksize = 7;

int main()
{
 char[N] a = "Hello ";
 int[N] b;
 b[0] = 15; b[1] = 10; b[2] = 6, b[3] = 0; b[4] = -11; b[5] = 1; b[6] = 0;

 char *ad;
 int *bd;
 const ulong csize = N*char.sizeof;
 const ulong isize = N*int.sizeof;

 writef("%s \n", a);

 cudaMalloc( cast(void**)&ad, csize );
 cudaMalloc( cast(void**)&bd, isize );
 cudaMemcpy( cast(void*)ad, cast(void*)a, csize, cudaMemcpyKind.cudaMemcpyHostToDevice );
 cudaMemcpy( cast(void*)bd, cast(void*)b, isize, cudaMemcpyKind.cudaMemcpyHostToDevice );

 hello(ad, bd, blocksize);
 cudaMemcpy(cast(void*)a, cast(void*)ad, csize, cudaMemcpyKind.cudaMemcpyDeviceToHost );
 cudaFree( ad );

 writef("%s \n", a);
 return 0;
}

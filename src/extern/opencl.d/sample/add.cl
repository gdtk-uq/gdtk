
#pragma OPENCL EXTENSION cl_khr_fp64 : enable


__kernel void set_array(
	__global float* src, float val)
{
	int gid = get_global_id(0);
	src[gid] = val;
}

__kernel void add_array(
	__global const float* src1,
	__global const float* src2,
	__global float* dst)
{
	int gid = get_global_id(0);
	dst[gid] = src1[gid] + src2[gid];
}

__kernel void set_array_double(
	__global double* src, double val)
{
	int gid = get_global_id(0);
	src[gid] = val;
}

__kernel void add_array_double(
	__global const double* src1,
	__global const double* src2,
	__global double* dst)
{
	int gid = get_global_id(0);
	dst[gid] = src1[gid] + src2[gid];
}


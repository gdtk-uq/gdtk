import std.array;
import std.stdio;
import std.conv;
import core.stdc.stdlib;
import core.stdc.stdio;
import core.memory;
import cl;

extern(C) void call_back(const char* e, const void* p, size_t cb, void* u)
{
        writefln("%s", e);
}

void main(string args[])
{
        const size_t VECTOR_SIZE = 64*1024*1024;

        const cl_uint num_entries = 10;
        cl_platform_id platform_ids[num_entries];
        cl_uint num_platforms;

        if(clGetPlatformIDs(
                num_entries, platform_ids.ptr, &num_platforms) != CL_SUCCESS)
        {
                exit(EXIT_FAILURE);
        }

        writefln("num_platforms: %d", num_platforms);

        for(cl_uint i = 0; (i < num_entries) && (i < num_platforms); ++i)
        {
                writefln("platform_ids[%d]: 0x%X", i, platform_ids[i]);
        }

        if(args.length < 2)
        {
                writefln("usage: %s kernel_file", args[0]);
        }

        auto kernel_file = File(args[1], "r");
        auto kernel_source = uninitializedArray!(char[])(kernel_file.size);
        kernel_source = kernel_file.rawRead(kernel_source);

        cl_device_id device_id = null;
        cl_uint ret_num_devices = 0;
        cl_int ret = clGetDeviceIDs(
                platform_ids[0], CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);

        cl_context context =
                clCreateContext(null, 1, &device_id, &call_back, null, &ret);
        cl_command_queue command_queue =
                clCreateCommandQueue(context, device_id, 0, &ret);

        size_t kernel_source_size = kernel_source.length;
        const char* p = kernel_source.ptr;
        cl_program program =
                clCreateProgramWithSource(context, 1, &p, &kernel_source_size, &ret);

        ret = clBuildProgram(program, 1, &device_id, null, null, null);

        size_t log_size = 0;
        clGetProgramBuildInfo(
                program, device_id, CL_PROGRAM_BUILD_LOG, 0, null, &log_size);

        char* log_buf = cast(char*)GC.malloc(log_size);
        clGetProgramBuildInfo(
                program, device_id, CL_PROGRAM_BUILD_LOG, log_size, &log_buf[0], &log_size);
        writefln("%s", to!string(log_buf));

        cl_kernel add_kernel = clCreateKernel(program, "add_array", &ret);
        cl_kernel set_kernel = clCreateKernel(program, "set_array", &ret);

        cl_mem memobj1 = clCreateBuffer(
                context, CL_MEM_READ_WRITE, VECTOR_SIZE * float.sizeof, null, &ret);
        cl_mem memobj2 = clCreateBuffer(
                context, CL_MEM_READ_WRITE, VECTOR_SIZE * float.sizeof, null, &ret);
        cl_mem memobj3 = clCreateBuffer(
                context, CL_MEM_READ_WRITE, VECTOR_SIZE * float.sizeof, null, &ret);

        float val1 = 1;
        float val2 = 2;

        cl_event events[3];

        size_t worksize[1] = [VECTOR_SIZE];
        size_t* localworksize = null;
        size_t* offset = null;
        
        ret = clSetKernelArg(set_kernel, 0, cl_mem.sizeof, cast(void*)&memobj1);
        ret = clSetKernelArg(set_kernel, 1, float.sizeof, cast(void*)(&val1));
        ret = clEnqueueNDRangeKernel(
                command_queue, set_kernel, 1,  offset, worksize.ptr, localworksize,
                0, null, &events[0]);

        ret = clSetKernelArg(set_kernel, 0, cl_mem.sizeof, cast(void*)&memobj2);
        ret = clSetKernelArg(set_kernel, 1, float.sizeof, cast(void*)(&val2));
        ret = clEnqueueNDRangeKernel(
                command_queue, set_kernel, 1,  offset, worksize.ptr, localworksize,
                0, null, &events[1]);

        float* resultFloats0 = cast(float*)GC.malloc(VECTOR_SIZE * float.sizeof);

        ret = clEnqueueReadBuffer(command_queue, memobj1, CL_TRUE, 0,
                VECTOR_SIZE * float.sizeof, resultFloats0, 2, &events[0], &events[2]);

        ret = clWaitForEvents(1, &events[0]);

        ret = clSetKernelArg(add_kernel, 0, cl_mem.sizeof, cast(void*)&memobj1);
        ret = clSetKernelArg(add_kernel, 1, cl_mem.sizeof, cast(void*)&memobj2);
        ret = clSetKernelArg(add_kernel, 2, cl_mem.sizeof, cast(void*)&memobj3);

        ret = clEnqueueNDRangeKernel(
                command_queue, add_kernel, 1,  null, worksize.ptr, null,
                2, &events[0], &events[2]);

        float* resultFloats3 = cast(float*)GC.malloc(VECTOR_SIZE * float.sizeof);

        ret = clEnqueueReadBuffer(command_queue, memobj3, CL_TRUE, 0,
                VECTOR_SIZE * float.sizeof, resultFloats3, 1, &events[0], &events[1]);
        ret = clWaitForEvents(1, &events[1]);

        for(int i = 0; i < 10; ++i)
        {
                writefln("%f", resultFloats3[i]);
        }

        clWaitForEvents(3, &events[0]);

        ret = clFlush(command_queue);
        ret = clFinish(command_queue);
        ret = clReleaseKernel(add_kernel);
        ret = clReleaseKernel(set_kernel);
        ret = clReleaseProgram(program);
        ret = clReleaseMemObject(memobj1);
        ret = clReleaseMemObject(memobj2);
        ret = clReleaseMemObject(memobj3);
        ret = clReleaseCommandQueue(command_queue);
        ret = clReleaseContext(context);
}


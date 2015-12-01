module cl;

import std.stdint;

alias byte cl_char;
alias ubyte cl_uchar;
alias short cl_short;
alias ushort cl_ushort;
alias int cl_int;
alias uint cl_uint;
alias long cl_long;
alias ulong cl_ulong;
alias ushort cl_half;
alias float cl_float;
alias double cl_double;

alias void* cl_platform_id;
alias void* cl_device_id;
alias void* cl_context;
alias void* cl_command_queue;
alias void* cl_mem;
alias void* cl_program;
alias void* cl_kernel;
alias void* cl_event;
alias void* cl_sampler;

alias cl_uint cl_bool;

alias cl_ulong cl_bitfield;
alias cl_bitfield cl_device_type;
alias cl_uint cl_platform_info;
alias cl_uint cl_device_info;
alias cl_bitfield cl_device_fp_config;
alias cl_uint cl_device_mem_cache_type;
alias cl_uint cl_device_local_mem_type;
alias cl_bitfield cl_device_exec_capabilities;
alias cl_bitfield cl_command_queue_properties;

alias intptr_t cl_context_properties;
alias cl_uint cl_context_info;
alias cl_uint cl_command_queue_info;
alias cl_uint cl_channel_order;
alias cl_uint cl_channel_type;
alias cl_bitfield cl_mem_flags;
alias cl_uint cl_mem_object_type;
alias cl_uint cl_mem_info;
alias cl_uint cl_image_info;
alias cl_uint cl_buffer_create_type;
alias cl_uint cl_addressing_mode;
alias cl_uint cl_filter_mode;
alias cl_uint cl_sampler_info;
alias cl_bitfield cl_map_flags;
alias cl_uint cl_program_info;
alias cl_uint cl_program_build_info;
alias cl_int cl_build_status;
alias cl_uint cl_kernel_info;
alias cl_uint cl_kernel_work_group_info;
alias cl_uint cl_event_info;
alias cl_uint cl_command_type;
alias cl_uint cl_profiling_info;

struct _cl_image_format
{
	cl_channel_order image_channel_order;
	cl_channel_type image;
}

alias _cl_image_format cl_image_format;

struct _cl_buffer_region
{
	size_t origin;
	size_t size;
}

alias _cl_buffer_region cl_buffer_region;

enum
{
	CL_SUCCESS = 0,
	CL_DEVICE_NOT_FOUND = -1,
	CL_DEVICE_NOT_AVAILABLE = -2,
	CL_COMPILER_NOT_AVAILABLE = -3,
	CL_MEM_OBJECT_ALLOCATION_FAILURE = -4,
	CL_OUT_OF_RESOURCES = -5,
	CL_OUT_OF_HOST_MEMORY = -6,
	CL_PROFILING_INFO_NOT_AVAILABLE = -7,
	CL_MEM_COPY_OVERLAP = -8,
	CL_IMAGE_FORMAT_MISMATCH = -9,
	CL_IMAGE_FORMAT_NOT_SUPPORTED = -10,
	CL_BUILD_PROGRAM_FAILURE = -11,
	CL_MAP_FAILURE = -12,
	CL_MISALIGNED_SUB_BUFFER_OFFSET = -13,
	CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST = -14,

	CL_INVALID_VALUE = -30,
	CL_INVALID_DEVICE_TYPE = -31,
	CL_INVALID_PLATFORM = -32,
	CL_INVALID_DEVICE = -33,
	CL_INVALID_CONTEXT = -34,
	CL_INVALID_QUEUE_PROPERTIES = -35,
	CL_INVALID_COMMAND_QUEUE = -36,
	CL_INVALID_HOST_PTR = -37,
	CL_INVALID_MEM_OBJECT = -38,
	CL_INVALID_IMAGE_FORMAT_DESCRIPTOR = -39,
	CL_INVALID_IMAGE_SIZE = -40,
	CL_INVALID_SAMPLER = -41,
	CL_INVALID_BINARY = -42,
	CL_INVALID_BUILD_OPTIONS = -43,
	CL_INVALID_PROGRAM = -44,
	CL_INVALID_PROGRAM_EXECUTABLE = -45,
	CL_INVALID_KERNEL_NAME = -46,
	CL_INVALID_KERNEL_DEFINITION = -47,
	CL_INVALID_KERNEL = -48,
	CL_INVALID_ARG_INDEX = -49,
	CL_INVALID_ARG_VALUE = -50,
	CL_INVALID_ARG_SIZE = -51,
	CL_INVALID_KERNEL_ARGS = -52,
	CL_INVALID_WORK_DIMENSION = -53,
	CL_INVALID_WORK_GROUP_SIZE = -54,
	CL_INVALID_WORK_ITEM_SIZE = -55,
	CL_INVALID_GLOBAL_OFFSET = -56,
	CL_INVALID_EVENT_WAIT_LIST = -57,
	CL_INVALID_EVENT = -58,
	CL_INVALID_OPERATION = -59,
	CL_INVALID_GL_OBJECT = -60,
	CL_INVALID_BUFFER_SIZE = -61,
	CL_INVALID_MIP_LEVEL = -62,
	CL_INVALID_GLOBAL_WORK_SIZE = -63,
	CL_INVALID_PROPERTY = -64
}

const int CL_VERSION_1_0 = 1;
const int CL_VERSION_1_1 = 1;

const int CL_FALSE = 0;
const int CL_TRUE = 1;

// cl_platform_info
const uint CL_PLATFORM_PROFILE = 0x0900;
const uint CL_PLATFORM_VERSION = 0x0901;
const uint CL_PLATFORM_NAME = 0x0902;
const uint CL_PLATFORM_VENDOR = 0x0903;
const uint CL_PLATFORM_EXTENSIONS = 0x0904;

// cl_device_type - bitfield
const uint CL_DEVICE_TYPE_DEFAULT = (1 << 0);
const uint CL_DEVICE_TYPE_CPU = (1 << 1);
const uint CL_DEVICE_TYPE_GPU = (1 << 2);
const uint CL_DEVICE_TYPE_ACCELERATOR = (1 << 3);
const uint CL_DEVICE_TYPE_ALL = 0xFFFFFFFF;

// cl_device_info
const uint CL_DEVICE_TYPE = 0x1000;
const uint CL_DEVICE_VENDOR_ID = 0x1001;
const uint CL_DEVICE_MAX_COMPUTE_UNITS = 0x1002;
const uint CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS = 0x1003;
const uint CL_DEVICE_MAX_WORK_GROUP_SIZE = 0x1004;
const uint CL_DEVICE_MAX_WORK_ITEM_SIZES = 0x1005;
const uint CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR = 0x1006;
const uint CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT = 0x1007;
const uint CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT = 0x1008;
const uint CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG = 0x1009;
const uint CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT = 0x100A;
const uint CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE = 0x100B;
const uint CL_DEVICE_MAX_CLOCK_FREQUENCY = 0x100C;
const uint CL_DEVICE_ADDRESS_BITS = 0x100D;
const uint CL_DEVICE_MAX_READ_IMAGE_ARGS = 0x100E;
const uint CL_DEVICE_MAX_WRITE_IMAGE_ARGS = 0x100F;
const uint CL_DEVICE_MAX_MEM_ALLOC_SIZE = 0x1010;
const uint CL_DEVICE_IMAGE2D_MAX_WIDTH = 0x1011;
const uint CL_DEVICE_IMAGE2D_MAX_HEIGHT = 0x1012;
const uint CL_DEVICE_IMAGE3D_MAX_WIDTH = 0x1013;
const uint CL_DEVICE_IMAGE3D_MAX_HEIGHT = 0x1014;
const uint CL_DEVICE_IMAGE3D_MAX_DEPTH = 0x1015;
const uint CL_DEVICE_IMAGE_SUPPORT = 0x1016;
const uint CL_DEVICE_MAX_PARAMETER_SIZE = 0x1017;
const uint CL_DEVICE_MAX_SAMPLERS = 0x1018;
const uint CL_DEVICE_MEM_BASE_ADDR_ALIGN = 0x1019;
const uint CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE = 0x101A;
const uint CL_DEVICE_SINGLE_FP_CONFIG = 0x101B;
const uint CL_DEVICE_GLOBAL_MEM_CACHE_TYPE = 0x101C;
const uint CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE = 0x101D;
const uint CL_DEVICE_GLOBAL_MEM_CACHE_SIZE = 0x101E;
const uint CL_DEVICE_GLOBAL_MEM_SIZE = 0x101F;
const uint CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE = 0x1020;
const uint CL_DEVICE_MAX_CONSTANT_ARGS = 0x1021;
const uint CL_DEVICE_LOCAL_MEM_TYPE = 0x1022;
const uint CL_DEVICE_LOCAL_MEM_SIZE = 0x1023;
const uint CL_DEVICE_ERROR_CORRECTION_SUPPORT = 0x1024;
const uint CL_DEVICE_PROFILING_TIMER_RESOLUTION = 0x1025;
const uint CL_DEVICE_ENDIAN_LITTLE = 0x1026;
const uint CL_DEVICE_AVAILABLE = 0x1027;
const uint CL_DEVICE_COMPILER_AVAILABLE = 0x1028;
const uint CL_DEVICE_EXECUTION_CAPABILITIES = 0x1029;
const uint CL_DEVICE_QUEUE_PROPERTIES = 0x102A;
const uint CL_DEVICE_NAME = 0x102B;
const uint CL_DEVICE_VENDOR = 0x102C;
const uint CL_DRIVER_VERSION = 0x102D;
const uint CL_DEVICE_PROFILE = 0x102E;
const uint CL_DEVICE_VERSION = 0x102F;
const uint CL_DEVICE_EXTENSIONS = 0x1030;
const uint CL_DEVICE_PLATFORM = 0x1031;
// 0x1032 reserved for CL_DEVICE_DOUBLE_FP_CONFIG
// 0x1033 reserved for CL_DEVICE_HALF_FP_CONFIG
const uint CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF = 0x1034;
const uint CL_DEVICE_HOST_UNIFIED_MEMORY = 0x1035;
const uint CL_DEVICE_NATIVE_VECTOR_WIDTH_CHAR = 0x1036;
const uint CL_DEVICE_NATIVE_VECTOR_WIDTH_SHORT = 0x1037;
const uint CL_DEVICE_NATIVE_VECTOR_WIDTH_INT = 0x1038;
const uint CL_DEVICE_NATIVE_VECTOR_WIDTH_LONG = 0x1039;
const uint CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT = 0x103A;
const uint CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE = 0x103B;
const uint CL_DEVICE_NATIVE_VECTOR_WIDTH_HALF = 0x103C;
const uint CL_DEVICE_OPENCL_C_VERSION = 0x103D;

// cl_device_fp_config - bitfield
const uint CL_FP_DENORM = (1 << 0);
const uint CL_FP_INF_NAN = (1 << 1);
const uint CL_FP_ROUND_TO_NEAREST = (1 << 2);
const uint CL_FP_ROUND_TO_ZERO = (1 << 3);
const uint CL_FP_ROUND_TO_INF = (1 << 4);
const uint CL_FP_FMA = (1 << 5);
const uint CL_FP_SOFT_FLOAT = (1 << 6);

// cl_device_mem_cache_type
const uint CL_NONE = 0x0;
const uint CL_READ_ONLY_CACHE = 0x1;
const uint CL_READ_WRITE_CACHE = 0x2;

// cl_device_local_mem_type
const uint CL_LOCAL = 0x1;
const uint CL_GLOBAL = 0x2;

// cl_device_exec_capabilities - bitfield
const uint CL_EXEC_KERNEL = (1 << 0);
const uint CL_EXEC_NATIVE_KERNEL = (1 << 1);

// cl_command_queue_properties - bitfield
const uint CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE = (1 << 0);
const uint CL_QUEUE_PROFILING_ENABLE = (1 << 1);

// cl_context_info
const uint CL_CONTEXT_REFERENCE_COUNT = 0x1080;
const uint CL_CONTEXT_DEVICES = 0x1081;
const uint CL_CONTEXT_PROPERTIES = 0x1082;
const uint CL_CONTEXT_NUM_DEVICES = 0x1083;

// cl_context_info + cl_context_properties
const uint CL_CONTEXT_PLATFORM = 0x1084;

// cl_command_queue_info
const uint CL_QUEUE_CONTEXT = 0x1090;
const uint CL_QUEUE_DEVICE = 0x1091;
const uint CL_QUEUE_REFERENCE_COUNT = 0x1092;
const uint CL_QUEUE_PROPERTIES = 0x1093;

// cl_mem_flags - bitfield
const uint CL_MEM_READ_WRITE = (1 << 0);
const uint CL_MEM_WRITE_ONLY = (1 << 1);
const uint CL_MEM_READ_ONLY = (1 << 2);
const uint CL_MEM_USE_HOST_PTR = (1 << 3);
const uint CL_MEM_ALLOC_HOST_PTR = (1 << 4);
const uint CL_MEM_COPY_HOST_PTR = (1 << 5);

// cl_channel_order
const uint CL_R = 0x10B0;
const uint CL_A = 0x10B1;
const uint CL_RG = 0x10B2;
const uint CL_RA = 0x10B3;
const uint CL_RGB = 0x10B4;
const uint CL_RGBA = 0x10B5;
const uint CL_BGRA = 0x10B6;
const uint CL_ARGB = 0x10B7;
const uint CL_INTENSITY = 0x10B8;
const uint CL_LUMINANCE = 0x10B9;
const uint CL_Rx = 0x10BA;
const uint CL_RGx = 0x10BB;
const uint CL_RGBx = 0x10BC;

// cl_channel_type
const uint CL_SNORM_INT8 = 0x10D0;
const uint CL_SNORM_INT16 = 0x10D1;
const uint CL_UNORM_INT8 = 0x10D2;
const uint CL_UNORM_INT16 = 0x10D3;
const uint CL_UNORM_SHORT_565 = 0x10D4;
const uint CL_UNORM_SHORT_555 = 0x10D5;
const uint CL_UNORM_INT_101010 = 0x10D6;
const uint CL_SIGNED_INT8 = 0x10D7;
const uint CL_SIGNED_INT16 = 0x10D8;
const uint CL_SIGNED_INT32 = 0x10D9;
const uint CL_UNSIGNED_INT8 = 0x10DA;
const uint CL_UNSIGNED_INT16 = 0x10DB;
const uint CL_UNSIGNED_INT32 = 0x10DC;
const uint CL_HALF_FLOAT = 0x10DD;
const uint CL_FLOAT = 0x10DE;

// cl_mem_object_type
const uint CL_MEM_OBJECT_BUFFER = 0x10F0;
const uint CL_MEM_OBJECT_IMAGE2D = 0x10F1;
const uint CL_MEM_OBJECT_IMAGE3D = 0x10F2;

// cl_mem_info
const uint CL_MEM_TYPE = 0x1100;
const uint CL_MEM_FLAGS = 0x1101;
const uint CL_MEM_SIZE = 0x1102;
const uint CL_MEM_HOST_PTR = 0x1103;
const uint CL_MEM_MAP_COUNT = 0x1104;
const uint CL_MEM_REFERENCE_COUNT = 0x1105;
const uint CL_MEM_CONTEXT = 0x1106;
const uint CL_MEM_ASSOCIATED_MEMOBJECT = 0x1107;
const uint CL_MEM_OFFSET = 0x1108;

// cl_image_info
const uint CL_IMAGE_FORMAT = 0x1110;
const uint CL_IMAGE_ELEMENT_SIZE = 0x1111;
const uint CL_IMAGE_ROW_PITCH = 0x1112;
const uint CL_IMAGE_SLICE_PITCH = 0x1113;
const uint CL_IMAGE_WIDTH = 0x1114;
const uint CL_IMAGE_HEIGHT = 0x1115;
const uint CL_IMAGE_DEPTH = 0x1116;

// cl_addressing_mode
const uint CL_ADDRESS_NONE = 0x1130;
const uint CL_ADDRESS_CLAMP_TO_EDGE = 0x1131;
const uint CL_ADDRESS_CLAMP = 0x1132;
const uint CL_ADDRESS_REPEAT = 0x1133;
const uint CL_ADDRESS_MIRRORED_REPEAT = 0x1134;

// cl_filter_mode
const uint CL_FILTER_NEAREST = 0x1140;
const uint CL_FILTER_LINEAR = 0x1141;

// cl_sampler_info
const uint CL_SAMPLER_REFERENCE_COUNT = 0x1150;
const uint CL_SAMPLER_CONTEXT = 0x1151;
const uint CL_SAMPLER_NORMALIZED_COORDS = 0x1152;
const uint CL_SAMPLER_ADDRESSING_MODE = 0x1153;
const uint CL_SAMPLER_FILTER_MODE = 0x1154;

// cl_map_flags - bitfield
const uint CL_MAP_READ = (1 << 0);
const uint CL_MAP_WRITE = (1 << 1);

// cl_program_info
const uint CL_PROGRAM_REFERENCE_COUNT = 0x1160;
const uint CL_PROGRAM_CONTEXT = 0x1161;
const uint CL_PROGRAM_NUM_DEVICES = 0x1162;
const uint CL_PROGRAM_DEVICES = 0x1163;
const uint CL_PROGRAM_SOURCE = 0x1164;
const uint CL_PROGRAM_BINARY_SIZES = 0x1165;
const uint CL_PROGRAM_BINARIES = 0x1166;

// cl_program_build_info
const uint CL_PROGRAM_BUILD_STATUS = 0x1181;
const uint CL_PROGRAM_BUILD_OPTIONS = 0x1182;
const uint CL_PROGRAM_BUILD_LOG = 0x1183;

// cl_build_status
const uint CL_BUILD_SUCCESS = 0;
const uint CL_BUILD_NONE = -1;
const uint CL_BUILD_ERROR = -2;
const uint CL_BUILD_IN_PROGRESS = -3;

// cl_kernel_info
const uint CL_KERNEL_FUNCTION_NAME = 0x1190;
const uint CL_KERNEL_NUM_ARGS = 0x1191;
const uint CL_KERNEL_REFERENCE_COUNT = 0x1192;
const uint CL_KERNEL_CONTEXT = 0x1193;
const uint CL_KERNEL_PROGRAM = 0x1194;

// cl_kernel_work_group_info
const uint CL_KERNEL_WORK_GROUP_SIZE = 0x11B0;
const uint CL_KERNEL_COMPILE_WORK_GROUP_SIZE = 0x11B1;
const uint CL_KERNEL_LOCAL_MEM_SIZE = 0x11B2;
const uint CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE = 0x11B3;
const uint CL_KERNEL_PRIVATE_MEM_SIZE = 0x11B4;

// cl_event_info
const uint CL_EVENT_COMMAND_QUEUE = 0x11D0;
const uint CL_EVENT_COMMAND_TYPE = 0x11D1;
const uint CL_EVENT_REFERENCE_COUNT = 0x11D2;
const uint CL_EVENT_COMMAND_EXECUTION_STATUS = 0x11D3;
const uint CL_EVENT_CONTEXT = 0x11D4;

// cl_command_type
const uint CL_COMMAND_NDRANGE_KERNEL = 0x11F0;
const uint CL_COMMAND_TASK = 0x11F1;
const uint CL_COMMAND_NATIVE_KERNEL = 0x11F2;
const uint CL_COMMAND_READ_BUFFER = 0x11F3;
const uint CL_COMMAND_WRITE_BUFFER = 0x11F4;
const uint CL_COMMAND_COPY_BUFFER = 0x11F5;
const uint CL_COMMAND_READ_IMAGE = 0x11F6;
const uint CL_COMMAND_WRITE_IMAGE = 0x11F7;
const uint CL_COMMAND_COPY_IMAGE = 0x11F8;
const uint CL_COMMAND_COPY_IMAGE_TO_BUFFER = 0x11F9;
const uint CL_COMMAND_COPY_BUFFER_TO_IMAGE = 0x11FA;
const uint CL_COMMAND_MAP_BUFFER = 0x11FB;
const uint CL_COMMAND_MAP_IMAGE = 0x11FC;
const uint CL_COMMAND_UNMAP_MEM_OBJECT = 0x11FD;
const uint CL_COMMAND_MARKER = 0x11FE;
const uint CL_COMMAND_ACQUIRE_GL_OBJECTS = 0x11FF;
const uint CL_COMMAND_RELEASE_GL_OBJECTS = 0x1200;
const uint CL_COMMAND_READ_BUFFER_RECT = 0x1201;
const uint CL_COMMAND_WRITE_BUFFER_RECT = 0x1202;
const uint CL_COMMAND_COPY_BUFFER_RECT = 0x1203;
const uint CL_COMMAND_USER = 0x1204;

// command execution status
const uint CL_COMPLETE = 0x0;
const uint CL_RUNNING = 0x1;
const uint CL_SUBMITTED = 0x2;
const uint CL_QUEUED = 0x3;

// cl_buffer_create_type
const uint CL_BUFFER_CREATE_TYPE_REGION = 0x1220;

// cl_profiling_info
const uint CL_PROFILING_COMMAND_QUEUED = 0x1280;
const uint CL_PROFILING_COMMAND_SUBMIT = 0x1281;
const uint CL_PROFILING_COMMAND_START = 0x1282;
const uint CL_PROFILING_COMMAND_END = 0x1283;

extern (System)
{

// Platform API
cl_int
clGetPlatformIDs(cl_uint, cl_platform_id *, cl_uint *);

cl_int 
clGetPlatformInfo(cl_platform_id, cl_platform_info, size_t, void *, size_t *);

// Device APIs 
cl_int
clGetDeviceIDs(cl_platform_id, cl_device_type, cl_uint, cl_device_id *, cl_uint *);

cl_int
clGetDeviceInfo(cl_device_id, cl_device_info, size_t, void *, size_t *);

// Context APIs 
cl_context
clCreateContext(const cl_context_properties *, cl_uint, const cl_device_id *,
	void function(const char *, const void *, size_t, void *), void *, cl_int *);

cl_context
clCreateContextFromType(const cl_context_properties *, cl_device_type,
	void function(const char *, const void *, size_t, void *), void *, cl_int *);

cl_int
clRetainContext(cl_context);

cl_int
clReleaseContext(cl_context);

cl_int
clGetContextInfo(cl_context, cl_context_info, size_t, void *, size_t *);

// Command Queue APIs
cl_command_queue
clCreateCommandQueue(cl_context, cl_device_id, cl_command_queue_properties, cl_int *);

cl_int
clRetainCommandQueue(cl_command_queue);

cl_int
clReleaseCommandQueue(cl_command_queue);

cl_int
clGetCommandQueueInfo(cl_command_queue, cl_command_queue_info, size_t, void *, size_t *);

deprecated
{
cl_int
clSetCommandQueueProperty(cl_command_queue, cl_command_queue_properties, cl_bool,
	cl_command_queue_properties *);
}

// Memory Object APIs
cl_mem
clCreateBuffer(cl_context, cl_mem_flags, size_t, void *, cl_int *);

cl_mem
clCreateSubBuffer(cl_mem, cl_mem_flags, cl_buffer_create_type, const void *,
	cl_int *);

cl_mem
clCreateImage2D(cl_context, cl_mem_flags, const cl_image_format *, size_t, size_t,
	size_t, void *, cl_int *);

cl_mem
clCreateImage3D(cl_context, cl_mem_flags, const cl_image_format *, size_t, size_t,
	size_t, size_t, size_t, void *, cl_int *);

cl_int
clRetainMemObject(cl_mem);

cl_int
clReleaseMemObject(cl_mem);

cl_int
clGetSupportedImageFormats(cl_context, cl_mem_flags, cl_mem_object_type, cl_uint,
	cl_image_format *, cl_uint *);

cl_int
clGetMemObjectInfo(cl_mem, cl_mem_info, size_t, void *, size_t *);

cl_int
clGetImageInfo(cl_mem, cl_image_info, size_t, void *, size_t *);

cl_int
clSetMemObjectDestructorCallback(  cl_mem, void function( cl_mem, void*), void *);  

// Sampler APIs 
cl_sampler
clCreateSampler(cl_context, cl_bool, cl_addressing_mode, cl_filter_mode, cl_int *);

cl_int
clRetainSampler(cl_sampler);

cl_int
clReleaseSampler(cl_sampler);

cl_int
clGetSamplerInfo(cl_sampler, cl_sampler_info, size_t, void *, size_t *);

// Program Object APIs 
cl_program
clCreateProgramWithSource(cl_context, cl_uint, const char **, const size_t *,
	cl_int *);

cl_program
clCreateProgramWithBinary(cl_context, cl_uint, const cl_device_id *, const size_t *,
	const ubyte **, cl_int *, cl_int *);

cl_int
clRetainProgram(cl_program);

cl_int
clReleaseProgram(cl_program);

cl_int
clBuildProgram(cl_program, cl_uint, const cl_device_id *, const char *,
	void function(cl_program, void *), void *);

cl_int
clUnloadCompiler();

cl_int
clGetProgramInfo(cl_program, cl_program_info, size_t, void *, size_t *);

cl_int
clGetProgramBuildInfo(cl_program, cl_device_id, cl_program_build_info, size_t,
	void *, size_t *);

// Kernel Object APIs
cl_kernel
clCreateKernel(cl_program, const char *, cl_int *);

cl_int
clCreateKernelsInProgram(cl_program, cl_uint, cl_kernel *, cl_uint *);

cl_int
clRetainKernel(cl_kernel);

cl_int
clReleaseKernel(cl_kernel);

cl_int
clSetKernelArg(cl_kernel, cl_uint, size_t, const void *);

cl_int
clGetKernelInfo(cl_kernel, cl_kernel_info, size_t, void *, size_t *);

cl_int
clGetKernelWorkGroupInfo(cl_kernel, cl_device_id, cl_kernel_work_group_info,
	size_t, void *, size_t *);

// Event Object APIs 
cl_int
clWaitForEvents(cl_uint, const cl_event *);

cl_int
clGetEventInfo(cl_event, cl_event_info, size_t, void *, size_t *);

cl_event
clCreateUserEvent(cl_context, cl_int *);               

cl_int
clRetainEvent(cl_event);

cl_int
clReleaseEvent(cl_event);

cl_int
clSetUserEventStatus(cl_event, cl_int);

cl_int
clSetEventCallback( cl_event, cl_int, void function(cl_event, cl_int, void *),
	void *);

// Profiling APIs 
cl_int
clGetEventProfilingInfo(cl_event, cl_profiling_info, size_t, void *, size_t *);

// Flush and Finish APIs
cl_int
clFlush(cl_command_queue);

cl_int
clFinish(cl_command_queue);

// Enqueued Commands APIs
cl_int
clEnqueueReadBuffer(cl_command_queue, cl_mem, cl_bool, size_t, size_t, void *,
	cl_uint, const cl_event *, cl_event *);

cl_int
clEnqueueReadBufferRect(cl_command_queue, cl_mem, cl_bool, const size_t *,
	const size_t *, const size_t *, size_t, size_t, size_t, size_t,
	void *, cl_uint, const cl_event *, cl_event *);

cl_int
clEnqueueWriteBuffer(cl_command_queue, cl_mem, cl_bool, size_t, size_t,
	const void *, cl_uint, const cl_event *, cl_event *);

cl_int
clEnqueueWriteBufferRect(cl_command_queue, cl_mem, cl_bool, const size_t *,
	const size_t *, const size_t *, size_t, size_t, size_t, size_t,
	const void *, cl_uint, const cl_event *, cl_event *);

cl_int
clEnqueueCopyBuffer(cl_command_queue, cl_mem, cl_mem, size_t, size_t, size_t,
	cl_uint, const cl_event *, cl_event *);

cl_int
clEnqueueCopyBufferRect(cl_command_queue, cl_mem, cl_mem, const size_t *,
	const size_t *, const size_t *, size_t, size_t, size_t, size_t, cl_uint,
	const cl_event *, cl_event *);

cl_int
clEnqueueReadImage(cl_command_queue, cl_mem, cl_bool, const size_t *,
	const size_t *, size_t, size_t, void *, cl_uint, const cl_event *, cl_event *);

cl_int
clEnqueueWriteImage(cl_command_queue, cl_mem, cl_bool, const size_t *,
	const size_t *, size_t, size_t, const void *, cl_uint, const cl_event *,
	cl_event *);

cl_int
clEnqueueCopyImage(cl_command_queue, cl_mem, cl_mem, const size_t *,
	const size_t *, const size_t *, cl_uint, const cl_event *, cl_event *);

cl_int
clEnqueueCopyImageToBuffer(cl_command_queue, cl_mem, cl_mem, const size_t *,
	const size_t *, size_t, cl_uint, const cl_event *, cl_event *);

cl_int
clEnqueueCopyBufferToImage(cl_command_queue, cl_mem, cl_mem, size_t,
	const size_t *, const size_t *, cl_uint, const cl_event *, cl_event *);

void *
clEnqueueMapBuffer(cl_command_queue, cl_mem, cl_bool, cl_map_flags, size_t,
	size_t, cl_uint, const cl_event *, cl_event *, cl_int *);

void *
clEnqueueMapImage(cl_command_queue, cl_mem, cl_bool, cl_map_flags, const size_t *,
	const size_t *, size_t *, size_t *, cl_uint, const cl_event *, cl_event *,
	cl_int *);

cl_int
clEnqueueUnmapMemObject(cl_command_queue, cl_mem, void *, cl_uint, const cl_event *,
	cl_event *);

cl_int
clEnqueueNDRangeKernel(cl_command_queue, cl_kernel, cl_uint, const size_t *,
	const size_t *, const size_t *, cl_uint, const cl_event *, cl_event *);

cl_int
clEnqueueTask(cl_command_queue, cl_kernel, cl_uint, const cl_event *, cl_event *);

cl_int
clEnqueueNativeKernel(cl_command_queue, void function(void *), void *, size_t,
	cl_uint, const cl_mem *, const void **, cl_uint, const cl_event *, cl_event *);

cl_int
clEnqueueMarker(cl_command_queue, cl_event *);

cl_int
clEnqueueWaitForEvents(cl_command_queue, cl_uint, const cl_event *);

cl_int
clEnqueueBarrier(cl_command_queue);

void * clGetExtensionFunctionAddress(const char *);

}


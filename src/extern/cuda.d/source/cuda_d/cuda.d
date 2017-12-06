/*
  D bindings for CUDA.
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/
module cuda_d.cuda;

import core.stdc.config;

extern (C):

/*DEVICE_BUILTIN*/

alias void cudaArray;
alias void* cudaIpcEventHandle_t;
alias void* cudaIpcMemHandle_t;
alias CUevent_st * cudaEvent_t;
alias enum cudaError cudaError_t;
alias enum cudaOutputMode cudaOutputMode_t;

alias cudaArray *  cudaArray_const_t;
alias cudaArray *  cudaArray_t;
//alias CUeglStreamConnection_st *  cudaEglStreamConnection;
//alias enumcudaError cudaError_t;
//alias CUevent_st *  cudaEvent_t;
//alias cudaGraphicsResource *  cudaGraphicsResource_t;
//alias cudaMipmappedArray *  cudaMipmappedArray_const_t;
//alias cudaMipmappedArray *  cudaMipmappedArray_t;
//alias enumcudaOutputMode cudaOutputMode_t;
//alias CUstream_st *  cudaStream_t;
//alias ulong  cudaSurfaceObject_t;
//alias ulong  cudaTextureObject_t;
//alias CUuuid_st  cudaUUID_t;

struct cudaDeviceProp{
int  ECCEnabled;
int  asyncEngineCount;
int  canMapHostMemory;
int  clockRate;
int  computeMode;
int  concurrentKernels;
int  concurrentManagedAccess;
int  deviceOverlap;
int  globalL1CacheSupported;
int  hostNativeAtomicSupported;
int  integrated;
int  isMultiGpuBoard;
int  kernelExecTimeoutEnabled;
int  l2CacheSize;
int  localL1CacheSupported;
int  major;
int  managedMemory;
int[3]  maxGridSize;
int  maxSurface1D;
int[2]  maxSurface1DLayered;
int[2]  maxSurface2D;
int[3]  maxSurface2DLayered;
int[3]  maxSurface3D;
int  maxSurfaceCubemap;
int[2]  maxSurfaceCubemapLayered;
int  maxTexture1D;
int[2]  maxTexture1DLayered;
int  maxTexture1DLinear;
int  maxTexture1DMipmap;
int[2]  maxTexture2D;
int[2]  maxTexture2DGather;
int[3]  maxTexture2DLayered;
int[3]  maxTexture2DLinear;
int[2]  maxTexture2DMipmap;
int[3]  maxTexture3D;
int[3]  maxTexture3DAlt;
int  maxTextureCubemap;
int[2]  maxTextureCubemapLayered;
int[3]  maxThreadsDim;
int  maxThreadsPerBlock;
int  maxThreadsPerMultiProcessor;
size_t  memPitch;
int  memoryBusWidth;
int  memoryClockRate;
int  minor;
int  multiGpuBoardGroupID;
int  multiProcessorCount;
char[256]  name;
int  pageableMemoryAccess;
int  pciBusID;
int  pciDeviceID;
int  pciDomainID;
int  regsPerBlock;
int  regsPerMultiprocessor;
size_t  sharedMemPerBlock;
size_t  sharedMemPerMultiprocessor;
int  singleToDoublePrecisionPerfRatio;
int  streamPrioritiesSupported;
size_t  surfaceAlignment;
int  tccDriver;
size_t  textureAlignment;
size_t  texturePitchAlignment;
size_t  totalConstMem;
size_t  totalGlobalMem;
int  unifiedAddressing;
int  warpSize;
}

enum cudaChannelFormatKind {
  cudaChannelFormatKindSigned = 0,
  cudaChannelFormatKindUnsigned = 1,
  cudaChannelFormatKindFloat = 2,
  cudaChannelFormatKindNone = 3
}
enum cudaComputeMode {
  cudaComputeModeDefault = 0,
  cudaComputeModeExclusive = 1,
  cudaComputeModeProhibited = 2,
  cudaComputeModeExclusiveProcess = 3
}
enum cudaError {
  cudaSuccess = 0,
  cudaErrorMissingConfiguration = 1,
  cudaErrorMemoryAllocation = 2,
  cudaErrorInitializationError = 3,
  cudaErrorLaunchFailure = 4,
  cudaErrorPriorLaunchFailure = 5,
  cudaErrorLaunchTimeout = 6,
  cudaErrorLaunchOutOfResources = 7,
  cudaErrorInvalidDeviceFunction = 8,
  cudaErrorInvalidConfiguration = 9,
  cudaErrorInvalidDevice = 10,
  cudaErrorInvalidValue = 11,
  cudaErrorInvalidPitchValue = 12,
  cudaErrorInvalidSymbol = 13,
  cudaErrorMapBufferObjectFailed = 14,
  cudaErrorUnmapBufferObjectFailed = 15,
  cudaErrorInvalidHostPointer = 16,
  cudaErrorInvalidDevicePointer = 17,
  cudaErrorInvalidTexture = 18,
  cudaErrorInvalidTextureBinding = 19,
  cudaErrorInvalidChannelDescriptor = 20,
  cudaErrorInvalidMemcpyDirection = 21,
  cudaErrorAddressOfConstant = 22,
  cudaErrorTextureFetchFailed = 23,
  cudaErrorTextureNotBound = 24,
  cudaErrorSynchronizationError = 25,
  cudaErrorInvalidFilterSetting = 26,
  cudaErrorInvalidNormSetting = 27,
  cudaErrorMixedDeviceExecution = 28,
  cudaErrorCudartUnloading = 29,
  cudaErrorUnknown = 30,
  cudaErrorNotYetImplemented = 31,
  cudaErrorMemoryValueTooLarge = 32,
  cudaErrorInvalidResourceHandle = 33,
  cudaErrorNotReady = 34,
  cudaErrorInsufficientDriver = 35,
  cudaErrorSetOnActiveProcess = 36,
  cudaErrorInvalidSurface = 37,
  cudaErrorNoDevice = 38,
  cudaErrorECCUncorrectable = 39,
  cudaErrorSharedObjectSymbolNotFound = 40,
  cudaErrorSharedObjectInitFailed = 41,
  cudaErrorUnsupportedLimit = 42,
  cudaErrorDuplicateVariableName = 43,
  cudaErrorDuplicateTextureName = 44,
  cudaErrorDuplicateSurfaceName = 45,
  cudaErrorDevicesUnavailable = 46,
  cudaErrorInvalidKernelImage = 47,
  cudaErrorNoKernelImageForDevice = 48,
  cudaErrorIncompatibleDriverContext = 49,
  cudaErrorPeerAccessAlreadyEnabled = 50,
  cudaErrorPeerAccessNotEnabled = 51,
  cudaErrorDeviceAlreadyInUse = 54,
  cudaErrorProfilerDisabled = 55,
  cudaErrorProfilerNotInitialized = 56,
  cudaErrorProfilerAlreadyStarted = 57,
  cudaErrorProfilerAlreadyStopped = 58,
  cudaErrorStartupFailure = 0x7f,
  cudaErrorApiFailureBase = 10000
}
enum cudaFuncCache {
  cudaFuncCachePreferNone = 0,
  cudaFuncCachePreferShared = 1,
  cudaFuncCachePreferL1 = 2
}
enum cudaGraphicsCubeFace {
  cudaGraphicsCubeFacePositiveX = 0x00,
  cudaGraphicsCubeFaceNegativeX = 0x01,
  cudaGraphicsCubeFacePositiveY = 0x02,
  cudaGraphicsCubeFaceNegativeY = 0x03,
  cudaGraphicsCubeFacePositiveZ = 0x04,
  cudaGraphicsCubeFaceNegativeZ = 0x05
}
enum cudaGraphicsMapFlags {
  cudaGraphicsMapFlagsNone = 0,
  cudaGraphicsMapFlagsReadOnly = 1,
  cudaGraphicsMapFlagsWriteDiscard = 2
}
enum cudaGraphicsRegisterFlags {
  cudaGraphicsRegisterFlagsNone = 0,
  cudaGraphicsRegisterFlagsReadOnly = 1,
  cudaGraphicsRegisterFlagsWriteDiscard = 2,
  cudaGraphicsRegisterFlagsSurfaceLoadStore = 4
}
enum cudaLimit {
  cudaLimitStackSize = 0x00,
  cudaLimitPrintfFifoSize = 0x01,
  cudaLimitMallocHeapSize = 0x02
}
enum cudaMemcpyKind {
  cudaMemcpyHostToHost = 0,
  cudaMemcpyHostToDevice = 1,
  cudaMemcpyDeviceToHost = 2,
  cudaMemcpyDeviceToDevice = 3,
  cudaMemcpyDefault = 4
}
enum cudaMemoryType {
  cudaMemoryTypeHost = 1,
  cudaMemoryTypeDevice = 2
}

enum cudaSharedMemConfig {
  cudaSharedMemBankSizeDefault = 0,
  cudaSharedMemBankSizeFourByte = 1,
  cudaSharedMemBankSizeEightByte = 2
}

enum cudaOutputMode{
  cudaKeyValuePair = 0x00,
  cudaCSV = 0x01,
}

/**
 * CUDA API versioning support
 */

enum __CUDA_API_VERSION = 7050;
/* CUDA_FORCE_API_VERSION */

extern (D) auto __CUDA_API_PTDS(T)(auto ref T api)
{
    return api;
}

extern (D) auto __CUDA_API_PTSZ(T)(auto ref T api)
{
    return api;
}

alias cuDeviceTotalMem = cuDeviceTotalMem_v2;
alias cuCtxCreate = cuCtxCreate_v2;
alias cuModuleGetGlobal = cuModuleGetGlobal_v2;
alias cuMemGetInfo = cuMemGetInfo_v2;
alias cuMemAlloc = cuMemAlloc_v2;
alias cuMemAllocPitch = cuMemAllocPitch_v2;
alias cuMemFree = cuMemFree_v2;
alias cuMemGetAddressRange = cuMemGetAddressRange_v2;
alias cuMemAllocHost = cuMemAllocHost_v2;
alias cuMemHostGetDevicePointer = cuMemHostGetDevicePointer_v2;
alias cuMemcpyHtoD = cuMemcpyHtoD_v2;
alias cuMemcpyDtoH = cuMemcpyDtoH_v2;
alias cuMemcpyDtoD = cuMemcpyDtoD_v2;
alias cuMemcpyDtoA = cuMemcpyDtoA_v2;
alias cuMemcpyAtoD = cuMemcpyAtoD_v2;
alias cuMemcpyHtoA = cuMemcpyHtoA_v2;
alias cuMemcpyAtoH = cuMemcpyAtoH_v2;
alias cuMemcpyAtoA = cuMemcpyAtoA_v2;
alias cuMemcpyHtoAAsync = cuMemcpyHtoAAsync_v2;
alias cuMemcpyAtoHAsync = cuMemcpyAtoHAsync_v2;
alias cuMemcpy2D = cuMemcpy2D_v2;
alias cuMemcpy2DUnaligned = cuMemcpy2DUnaligned_v2;
alias cuMemcpy3D = cuMemcpy3D_v2;
alias cuMemcpyHtoDAsync = cuMemcpyHtoDAsync_v2;
alias cuMemcpyDtoHAsync = cuMemcpyDtoHAsync_v2;
alias cuMemcpyDtoDAsync = cuMemcpyDtoDAsync_v2;
alias cuMemcpy2DAsync = cuMemcpy2DAsync_v2;
alias cuMemcpy3DAsync = cuMemcpy3DAsync_v2;
alias cuMemsetD8 = cuMemsetD8_v2;
alias cuMemsetD16 = cuMemsetD16_v2;
alias cuMemsetD32 = cuMemsetD32_v2;
alias cuMemsetD2D8 = cuMemsetD2D8_v2;
alias cuMemsetD2D16 = cuMemsetD2D16_v2;
alias cuMemsetD2D32 =  cuMemsetD2D32_v2;
alias cuArrayCreate = cuArrayCreate_v2;
alias cuArrayGetDescriptor = cuArrayGetDescriptor_v2;
alias cuArray3DCreate = cuArray3DCreate_v2;
alias cuArray3DGetDescriptor = cuArray3DGetDescriptor_v2;
alias cuTexRefSetAddress = cuTexRefSetAddress_v2;
alias cuTexRefGetAddress = cuTexRefGetAddress_v2;
alias cuGraphicsResourceGetMappedPointer = cuGraphicsResourceGetMappedPointer_v2;
/* __CUDA_API_VERSION_INTERNAL || __CUDA_API_VERSION >= 3020 */
alias cuCtxDestroy = cuCtxDestroy_v2;
alias cuCtxPopCurrent = cuCtxPopCurrent_v2;
alias cuCtxPushCurrent = cuCtxPushCurrent_v2;
alias cuStreamDestroy = cuStreamDestroy_v2;
alias cuEventDestroy = cuEventDestroy_v2;
/* __CUDA_API_VERSION_INTERNAL || __CUDA_API_VERSION >= 4000 */
alias cuTexRefSetAddress2D = cuTexRefSetAddress2D_v3;
/* __CUDA_API_VERSION_INTERNAL || __CUDA_API_VERSION >= 4010 */
alias cuLinkCreate = cuLinkCreate_v2;
alias cuLinkAddData = cuLinkAddData_v2;
alias cuLinkAddFile = cuLinkAddFile_v2;
/* __CUDA_API_VERSION_INTERNAL || __CUDA_API_VERSION >= 6050 */
alias cuMemHostRegister = cuMemHostRegister_v2;
alias cuGraphicsResourceSetMapFlags = cuGraphicsResourceSetMapFlags_v2;

enum CUDA_VERSION = 7050;

/**
 * CUDA device pointer
 * CUdeviceptr is defined as an unsigned integer type whose size matches the size of a pointer on the target platform.
 */

alias CUdeviceptr = ulong;

/* __CUDA_API_VERSION >= 3020 */

alias CUdevice = int; /**< CUDA device */
struct CUctx_st;
alias CUcontext = CUctx_st*; /**< CUDA context */
struct CUmod_st;
alias CUmodule = CUmod_st*; /**< CUDA module */
struct CUfunc_st;
alias CUfunction = CUfunc_st*; /**< CUDA function */
struct CUarray_st;
alias CUarray = CUarray_st*; /**< CUDA array */
struct CUmipmappedArray_st;
alias CUmipmappedArray = CUmipmappedArray_st*; /**< CUDA mipmapped array */
struct CUtexref_st;
alias CUtexref = CUtexref_st*; /**< CUDA texture reference */
struct CUsurfref_st;
alias CUsurfref = CUsurfref_st*; /**< CUDA surface reference */
struct CUevent_st;
alias CUevent = CUevent_st*; /**< CUDA event */
struct CUstream_st;
alias CUstream = CUstream_st*; /**< CUDA stream */
struct CUgraphicsResource_st;
alias CUgraphicsResource = CUgraphicsResource_st*; /**< CUDA graphics interop resource */
alias CUtexObject = ulong; /**< An opaque value that represents a CUDA texture object */
alias CUsurfObject = ulong; /**< An opaque value that represents a CUDA surface object */

struct CUuuid_st
{
    /**< CUDA definition of UUID */
    char[16] bytes;
}

alias CUuuid = CUuuid_st;

/**
 * CUDA IPC handle size
 */
enum CU_IPC_HANDLE_SIZE = 64;

/**
 * CUDA IPC event handle
 */
struct CUipcEventHandle_st
{
    char[CU_IPC_HANDLE_SIZE] reserved;
}

alias CUipcEventHandle = CUipcEventHandle_st;

/**
 * CUDA IPC mem handle
 */
struct CUipcMemHandle_st
{
    char[CU_IPC_HANDLE_SIZE] reserved;
}

alias CUipcMemHandle = CUipcMemHandle_st;

/**
 * CUDA Ipc Mem Flags
 */
enum CUipcMem_flags_enum
{
    CU_IPC_MEM_LAZY_ENABLE_PEER_ACCESS = 1 /**< Automatically enable peer access between remote devices as needed */
}

alias CUipcMem_flags = CUipcMem_flags_enum;

/**
 * CUDA Mem Attach Flags
 */
enum CUmemAttach_flags_enum
{
    CU_MEM_ATTACH_GLOBAL = 1, /**< Memory can be accessed by any stream on any device */
    CU_MEM_ATTACH_HOST = 2, /**< Memory cannot be accessed by any stream on any device */
    CU_MEM_ATTACH_SINGLE = 4 /**< Memory can only be accessed by a single stream on the associated device */
}

alias CUmemAttach_flags = CUmemAttach_flags_enum;

/**
 * Context creation flags
 */
enum CUctx_flags_enum
{
    CU_CTX_SCHED_AUTO = 0, /**< Automatic scheduling */
    CU_CTX_SCHED_SPIN = 1, /**< Set spin as default scheduling */
    CU_CTX_SCHED_YIELD = 2, /**< Set yield as default scheduling */
    CU_CTX_SCHED_BLOCKING_SYNC = 4, /**< Set blocking synchronization as default scheduling */
    CU_CTX_BLOCKING_SYNC = 4, /**< Set blocking synchronization as default scheduling
      *  \deprecated This flag was deprecated as of CUDA 4.0
      *  and was replaced with ::CU_CTX_SCHED_BLOCKING_SYNC. */
    CU_CTX_SCHED_MASK = 7,
    CU_CTX_MAP_HOST = 8, /**< Support mapped pinned allocations */
    CU_CTX_LMEM_RESIZE_TO_MAX = 16, /**< Keep local memory allocation after launch */
    CU_CTX_FLAGS_MASK = 31
}

alias CUctx_flags = CUctx_flags_enum;

/**
 * Stream creation flags
 */
enum CUstream_flags_enum
{
    CU_STREAM_DEFAULT = 0, /**< Default stream flag */
    CU_STREAM_NON_BLOCKING = 1 /**< Stream does not synchronize with stream 0 (the NULL stream) */
}

alias CUstream_flags = CUstream_flags_enum;

/**
 * Legacy stream handle
 *
 * Stream handle that can be passed as a CUstream to use an implicit stream
 * with legacy synchronization behavior.
 *
 * See details of the \link_sync_behavior
 */
enum CU_STREAM_LEGACY = cast(CUstream) 0x1;

/**
 * Per-thread stream handle
 *
 * Stream handle that can be passed as a CUstream to use an implicit stream
 * with per-thread synchronization behavior.
 *
 * See details of the \link_sync_behavior
 */
enum CU_STREAM_PER_THREAD = cast(CUstream) 0x2;

/**
 * Event creation flags
 */
enum CUevent_flags_enum
{
    CU_EVENT_DEFAULT = 0, /**< Default event flag */
    CU_EVENT_BLOCKING_SYNC = 1, /**< Event uses blocking synchronization */
    CU_EVENT_DISABLE_TIMING = 2, /**< Event will not record timing data */
    CU_EVENT_INTERPROCESS = 4 /**< Event is suitable for interprocess use. CU_EVENT_DISABLE_TIMING must be set */
}

alias CUevent_flags = CUevent_flags_enum;

/**
 * Occupancy calculator flag
 */
enum CUoccupancy_flags_enum
{
    CU_OCCUPANCY_DEFAULT = 0, /**< Default behavior */
    CU_OCCUPANCY_DISABLE_CACHING_OVERRIDE = 1 /**< Assume global caching is enabled and cannot be automatically turned off */
}

alias CUoccupancy_flags = CUoccupancy_flags_enum;

/**
 * Array formats
 */
enum CUarray_format_enum
{
    CU_AD_FORMAT_UNSIGNED_INT8 = 1, /**< Unsigned 8-bit integers */
    CU_AD_FORMAT_UNSIGNED_INT16 = 2, /**< Unsigned 16-bit integers */
    CU_AD_FORMAT_UNSIGNED_INT32 = 3, /**< Unsigned 32-bit integers */
    CU_AD_FORMAT_SIGNED_INT8 = 8, /**< Signed 8-bit integers */
    CU_AD_FORMAT_SIGNED_INT16 = 9, /**< Signed 16-bit integers */
    CU_AD_FORMAT_SIGNED_INT32 = 10, /**< Signed 32-bit integers */
    CU_AD_FORMAT_HALF = 16, /**< 16-bit floating point */
    CU_AD_FORMAT_FLOAT = 32 /**< 32-bit floating point */
}

alias CUarray_format = CUarray_format_enum;

/**
 * Texture reference addressing modes
 */
enum CUaddress_mode_enum
{
    CU_TR_ADDRESS_MODE_WRAP = 0, /**< Wrapping address mode */
    CU_TR_ADDRESS_MODE_CLAMP = 1, /**< Clamp to edge address mode */
    CU_TR_ADDRESS_MODE_MIRROR = 2, /**< Mirror address mode */
    CU_TR_ADDRESS_MODE_BORDER = 3 /**< Border address mode */
}

alias CUaddress_mode = CUaddress_mode_enum;

/**
 * Texture reference filtering modes
 */
enum CUfilter_mode_enum
{
    CU_TR_FILTER_MODE_POINT = 0, /**< Point filter mode */
    CU_TR_FILTER_MODE_LINEAR = 1 /**< Linear filter mode */
}

alias CUfilter_mode = CUfilter_mode_enum;

/**
 * Device properties
 */
enum CUdevice_attribute_enum
{
    CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK = 1, /**< Maximum number of threads per block */
    CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X = 2, /**< Maximum block dimension X */
    CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Y = 3, /**< Maximum block dimension Y */
    CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Z = 4, /**< Maximum block dimension Z */
    CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_X = 5, /**< Maximum grid dimension X */
    CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Y = 6, /**< Maximum grid dimension Y */
    CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Z = 7, /**< Maximum grid dimension Z */
    CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK = 8, /**< Maximum shared memory available per block in bytes */
    CU_DEVICE_ATTRIBUTE_SHARED_MEMORY_PER_BLOCK = 8, /**< Deprecated, use CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK */
    CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY = 9, /**< Memory available on device for __constant__ variables in a CUDA C kernel in bytes */
    CU_DEVICE_ATTRIBUTE_WARP_SIZE = 10, /**< Warp size in threads */
    CU_DEVICE_ATTRIBUTE_MAX_PITCH = 11, /**< Maximum pitch in bytes allowed by memory copies */
    CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK = 12, /**< Maximum number of 32-bit registers available per block */
    CU_DEVICE_ATTRIBUTE_REGISTERS_PER_BLOCK = 12, /**< Deprecated, use CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK */
    CU_DEVICE_ATTRIBUTE_CLOCK_RATE = 13, /**< Typical clock frequency in kilohertz */
    CU_DEVICE_ATTRIBUTE_TEXTURE_ALIGNMENT = 14, /**< Alignment requirement for textures */
    CU_DEVICE_ATTRIBUTE_GPU_OVERLAP = 15, /**< Device can possibly copy memory and execute a kernel concurrently. Deprecated. Use instead CU_DEVICE_ATTRIBUTE_ASYNC_ENGINE_COUNT. */
    CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT = 16, /**< Number of multiprocessors on device */
    CU_DEVICE_ATTRIBUTE_KERNEL_EXEC_TIMEOUT = 17, /**< Specifies whether there is a run time limit on kernels */
    CU_DEVICE_ATTRIBUTE_INTEGRATED = 18, /**< Device is integrated with host memory */
    CU_DEVICE_ATTRIBUTE_CAN_MAP_HOST_MEMORY = 19, /**< Device can map host memory into CUDA address space */
    CU_DEVICE_ATTRIBUTE_COMPUTE_MODE = 20, /**< Compute mode (See ::CUcomputemode for details) */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE1D_WIDTH = 21, /**< Maximum 1D texture width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_WIDTH = 22, /**< Maximum 2D texture width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_HEIGHT = 23, /**< Maximum 2D texture height */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_WIDTH = 24, /**< Maximum 3D texture width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_HEIGHT = 25, /**< Maximum 3D texture height */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_DEPTH = 26, /**< Maximum 3D texture depth */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_LAYERED_WIDTH = 27, /**< Maximum 2D layered texture width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_LAYERED_HEIGHT = 28, /**< Maximum 2D layered texture height */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_LAYERED_LAYERS = 29, /**< Maximum layers in a 2D layered texture */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_ARRAY_WIDTH = 27, /**< Deprecated, use CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_LAYERED_WIDTH */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_ARRAY_HEIGHT = 28, /**< Deprecated, use CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_LAYERED_HEIGHT */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_ARRAY_NUMSLICES = 29, /**< Deprecated, use CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_LAYERED_LAYERS */
    CU_DEVICE_ATTRIBUTE_SURFACE_ALIGNMENT = 30, /**< Alignment requirement for surfaces */
    CU_DEVICE_ATTRIBUTE_CONCURRENT_KERNELS = 31, /**< Device can possibly execute multiple kernels concurrently */
    CU_DEVICE_ATTRIBUTE_ECC_ENABLED = 32, /**< Device has ECC support enabled */
    CU_DEVICE_ATTRIBUTE_PCI_BUS_ID = 33, /**< PCI bus ID of the device */
    CU_DEVICE_ATTRIBUTE_PCI_DEVICE_ID = 34, /**< PCI device ID of the device */
    CU_DEVICE_ATTRIBUTE_TCC_DRIVER = 35, /**< Device is using TCC driver model */
    CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE = 36, /**< Peak memory clock frequency in kilohertz */
    CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH = 37, /**< Global memory bus width in bits */
    CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE = 38, /**< Size of L2 cache in bytes */
    CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_MULTIPROCESSOR = 39, /**< Maximum resident threads per multiprocessor */
    CU_DEVICE_ATTRIBUTE_ASYNC_ENGINE_COUNT = 40, /**< Number of asynchronous engines */
    CU_DEVICE_ATTRIBUTE_UNIFIED_ADDRESSING = 41, /**< Device shares a unified address space with the host */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE1D_LAYERED_WIDTH = 42, /**< Maximum 1D layered texture width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE1D_LAYERED_LAYERS = 43, /**< Maximum layers in a 1D layered texture */
    CU_DEVICE_ATTRIBUTE_CAN_TEX2D_GATHER = 44, /**< Deprecated, do not use. */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_GATHER_WIDTH = 45, /**< Maximum 2D texture width if CUDA_ARRAY3D_TEXTURE_GATHER is set */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_GATHER_HEIGHT = 46, /**< Maximum 2D texture height if CUDA_ARRAY3D_TEXTURE_GATHER is set */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_WIDTH_ALTERNATE = 47, /**< Alternate maximum 3D texture width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_HEIGHT_ALTERNATE = 48, /**< Alternate maximum 3D texture height */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_DEPTH_ALTERNATE = 49, /**< Alternate maximum 3D texture depth */
    CU_DEVICE_ATTRIBUTE_PCI_DOMAIN_ID = 50, /**< PCI domain ID of the device */
    CU_DEVICE_ATTRIBUTE_TEXTURE_PITCH_ALIGNMENT = 51, /**< Pitch alignment requirement for textures */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURECUBEMAP_WIDTH = 52, /**< Maximum cubemap texture width/height */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURECUBEMAP_LAYERED_WIDTH = 53, /**< Maximum cubemap layered texture width/height */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURECUBEMAP_LAYERED_LAYERS = 54, /**< Maximum layers in a cubemap layered texture */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACE1D_WIDTH = 55, /**< Maximum 1D surface width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACE2D_WIDTH = 56, /**< Maximum 2D surface width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACE2D_HEIGHT = 57, /**< Maximum 2D surface height */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACE3D_WIDTH = 58, /**< Maximum 3D surface width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACE3D_HEIGHT = 59, /**< Maximum 3D surface height */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACE3D_DEPTH = 60, /**< Maximum 3D surface depth */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACE1D_LAYERED_WIDTH = 61, /**< Maximum 1D layered surface width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACE1D_LAYERED_LAYERS = 62, /**< Maximum layers in a 1D layered surface */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACE2D_LAYERED_WIDTH = 63, /**< Maximum 2D layered surface width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACE2D_LAYERED_HEIGHT = 64, /**< Maximum 2D layered surface height */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACE2D_LAYERED_LAYERS = 65, /**< Maximum layers in a 2D layered surface */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACECUBEMAP_WIDTH = 66, /**< Maximum cubemap surface width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACECUBEMAP_LAYERED_WIDTH = 67, /**< Maximum cubemap layered surface width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_SURFACECUBEMAP_LAYERED_LAYERS = 68, /**< Maximum layers in a cubemap layered surface */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE1D_LINEAR_WIDTH = 69, /**< Maximum 1D linear texture width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_LINEAR_WIDTH = 70, /**< Maximum 2D linear texture width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_LINEAR_HEIGHT = 71, /**< Maximum 2D linear texture height */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_LINEAR_PITCH = 72, /**< Maximum 2D linear texture pitch in bytes */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_MIPMAPPED_WIDTH = 73, /**< Maximum mipmapped 2D texture width */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_MIPMAPPED_HEIGHT = 74, /**< Maximum mipmapped 2D texture height */
    CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR = 75, /**< Major compute capability version number */
    CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR = 76, /**< Minor compute capability version number */
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE1D_MIPMAPPED_WIDTH = 77, /**< Maximum mipmapped 1D texture width */
    CU_DEVICE_ATTRIBUTE_STREAM_PRIORITIES_SUPPORTED = 78, /**< Device supports stream priorities */
    CU_DEVICE_ATTRIBUTE_GLOBAL_L1_CACHE_SUPPORTED = 79, /**< Device supports caching globals in L1 */
    CU_DEVICE_ATTRIBUTE_LOCAL_L1_CACHE_SUPPORTED = 80, /**< Device supports caching locals in L1 */
    CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_MULTIPROCESSOR = 81, /**< Maximum shared memory available per multiprocessor in bytes */
    CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_MULTIPROCESSOR = 82, /**< Maximum number of 32-bit registers available per multiprocessor */
    CU_DEVICE_ATTRIBUTE_MANAGED_MEMORY = 83, /**< Device can allocate managed memory on this system */
    CU_DEVICE_ATTRIBUTE_MULTI_GPU_BOARD = 84, /**< Device is on a multi-GPU board */
    CU_DEVICE_ATTRIBUTE_MULTI_GPU_BOARD_GROUP_ID = 85, /**< Unique id for a group of devices on the same multi-GPU board */
    CU_DEVICE_ATTRIBUTE_MAX = 86
}

alias CUdevice_attribute = CUdevice_attribute_enum;

/**
 * Legacy device properties
 */
struct CUdevprop_st
{
    int maxThreadsPerBlock; /**< Maximum number of threads per block */
    int[3] maxThreadsDim; /**< Maximum size of each dimension of a block */
    int[3] maxGridSize; /**< Maximum size of each dimension of a grid */
    int sharedMemPerBlock; /**< Shared memory available per block in bytes */
    int totalConstantMemory; /**< Constant memory available on device in bytes */
    int SIMDWidth; /**< Warp size in threads */
    int memPitch; /**< Maximum pitch in bytes allowed by memory copies */
    int regsPerBlock; /**< 32-bit registers available per block */
    int clockRate; /**< Clock frequency in kilohertz */
    int textureAlign; /**< Alignment requirement for textures */
}

alias CUdevprop = CUdevprop_st;

/**
 * Pointer information
 */
enum CUpointer_attribute_enum
{
    CU_POINTER_ATTRIBUTE_CONTEXT = 1, /**< The ::CUcontext on which a pointer was allocated or registered */
    CU_POINTER_ATTRIBUTE_MEMORY_TYPE = 2, /**< The ::CUmemorytype describing the physical location of a pointer */
    CU_POINTER_ATTRIBUTE_DEVICE_POINTER = 3, /**< The address at which a pointer's memory may be accessed on the device */
    CU_POINTER_ATTRIBUTE_HOST_POINTER = 4, /**< The address at which a pointer's memory may be accessed on the host */
    CU_POINTER_ATTRIBUTE_P2P_TOKENS = 5, /**< A pair of tokens for use with the nv-p2p.h Linux kernel interface */
    CU_POINTER_ATTRIBUTE_SYNC_MEMOPS = 6, /**< Synchronize every synchronous memory operation initiated on this region */
    CU_POINTER_ATTRIBUTE_BUFFER_ID = 7, /**< A process-wide unique ID for an allocated memory region*/
    CU_POINTER_ATTRIBUTE_IS_MANAGED = 8 /**< Indicates if the pointer points to managed memory */
}

alias CUpointer_attribute = CUpointer_attribute_enum;

/**
 * Function properties
 */
enum CUfunction_attribute_enum
{
    /**
     * The maximum number of threads per block, beyond which a launch of the
     * function would fail. This number depends on both the function and the
     * device on which the function is currently loaded.
     */
    CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK = 0,

    /**
     * The size in bytes of statically-allocated shared memory required by
     * this function. This does not include dynamically-allocated shared
     * memory requested by the user at runtime.
     */
    CU_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES = 1,

    /**
     * The size in bytes of user-allocated constant memory required by this
     * function.
     */
    CU_FUNC_ATTRIBUTE_CONST_SIZE_BYTES = 2,

    /**
     * The size in bytes of local memory used by each thread of this function.
     */
    CU_FUNC_ATTRIBUTE_LOCAL_SIZE_BYTES = 3,

    /**
     * The number of registers used by each thread of this function.
     */
    CU_FUNC_ATTRIBUTE_NUM_REGS = 4,

    /**
     * The PTX virtual architecture version for which the function was
     * compiled. This value is the major PTX version * 10 + the minor PTX
     * version, so a PTX version 1.3 function would return the value 13.
     * Note that this may return the undefined value of 0 for cubins
     * compiled prior to CUDA 3.0.
     */
    CU_FUNC_ATTRIBUTE_PTX_VERSION = 5,

    /**
     * The binary architecture version for which the function was compiled.
     * This value is the major binary version * 10 + the minor binary version,
     * so a binary version 1.3 function would return the value 13. Note that
     * this will return a value of 10 for legacy cubins that do not have a
     * properly-encoded binary architecture version.
     */
    CU_FUNC_ATTRIBUTE_BINARY_VERSION = 6,

    /**
     * The attribute to indicate whether the function has been compiled with
     * user specified option "-Xptxas --dlcm=ca" set .
     */
    CU_FUNC_ATTRIBUTE_CACHE_MODE_CA = 7,

    CU_FUNC_ATTRIBUTE_MAX = 8
}

alias CUfunction_attribute = CUfunction_attribute_enum;

/**
 * Function cache configurations
 */
enum CUfunc_cache_enum
{
    CU_FUNC_CACHE_PREFER_NONE = 0, /**< no preference for shared memory or L1 (default) */
    CU_FUNC_CACHE_PREFER_SHARED = 1, /**< prefer larger shared memory and smaller L1 cache */
    CU_FUNC_CACHE_PREFER_L1 = 2, /**< prefer larger L1 cache and smaller shared memory */
    CU_FUNC_CACHE_PREFER_EQUAL = 3 /**< prefer equal sized L1 cache and shared memory */
}

alias CUfunc_cache = CUfunc_cache_enum;

/**
 * Shared memory configurations
 */
enum CUsharedconfig_enum
{
    CU_SHARED_MEM_CONFIG_DEFAULT_BANK_SIZE = 0, /**< set default shared memory bank size */
    CU_SHARED_MEM_CONFIG_FOUR_BYTE_BANK_SIZE = 1, /**< set shared memory bank width to four bytes */
    CU_SHARED_MEM_CONFIG_EIGHT_BYTE_BANK_SIZE = 2 /**< set shared memory bank width to eight bytes */
}

alias CUsharedconfig = CUsharedconfig_enum;

/**
 * Memory types
 */
enum CUmemorytype_enum
{
    CU_MEMORYTYPE_HOST = 1, /**< Host memory */
    CU_MEMORYTYPE_DEVICE = 2, /**< Device memory */
    CU_MEMORYTYPE_ARRAY = 3, /**< Array memory */
    CU_MEMORYTYPE_UNIFIED = 4 /**< Unified device or host memory */
}

alias CUmemorytype = CUmemorytype_enum;

/**
 * Compute Modes
 */
enum CUcomputemode_enum
{
    CU_COMPUTEMODE_DEFAULT = 0, /**< Default compute mode (Multiple contexts allowed per device) */
    CU_COMPUTEMODE_EXCLUSIVE = 1, /**< Compute-exclusive-thread mode (Only one context used by a single thread can be present on this device at a time) */
    CU_COMPUTEMODE_PROHIBITED = 2, /**< Compute-prohibited mode (No contexts can be created on this device at this time) */
    CU_COMPUTEMODE_EXCLUSIVE_PROCESS = 3 /**< Compute-exclusive-process mode (Only one context used by a single process can be present on this device at a time) */
}

alias CUcomputemode = CUcomputemode_enum;

/**
 * Online compiler and linker options
 */
enum CUjit_option_enum
{
    /**
     * Max number of registers that a thread may use.\n
     * Option type: unsigned int\n
     * Applies to: compiler only
     */
    CU_JIT_MAX_REGISTERS = 0,

    /**
     * IN: Specifies minimum number of threads per block to target compilation
     * for\n
     * OUT: Returns the number of threads the compiler actually targeted.
     * This restricts the resource utilization fo the compiler (e.g. max
     * registers) such that a block with the given number of threads should be
     * able to launch based on register limitations. Note, this option does not
     * currently take into account any other resource limitations, such as
     * shared memory utilization.\n
     * Cannot be combined with ::CU_JIT_TARGET.\n
     * Option type: unsigned int\n
     * Applies to: compiler only
     */
    CU_JIT_THREADS_PER_BLOCK = 1,

    /**
     * Overwrites the option value with the total wall clock time, in
     * milliseconds, spent in the compiler and linker\n
     * Option type: float\n
     * Applies to: compiler and linker
     */
    CU_JIT_WALL_TIME = 2,

    /**
     * Pointer to a buffer in which to print any log messages
     * that are informational in nature (the buffer size is specified via
     * option ::CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES)\n
     * Option type: char *\n
     * Applies to: compiler and linker
     */
    CU_JIT_INFO_LOG_BUFFER = 3,

    /**
     * IN: Log buffer size in bytes.  Log messages will be capped at this size
     * (including null terminator)\n
     * OUT: Amount of log buffer filled with messages\n
     * Option type: unsigned int\n
     * Applies to: compiler and linker
     */
    CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES = 4,

    /**
     * Pointer to a buffer in which to print any log messages that
     * reflect errors (the buffer size is specified via option
     * ::CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES)\n
     * Option type: char *\n
     * Applies to: compiler and linker
     */
    CU_JIT_ERROR_LOG_BUFFER = 5,

    /**
     * IN: Log buffer size in bytes.  Log messages will be capped at this size
     * (including null terminator)\n
     * OUT: Amount of log buffer filled with messages\n
     * Option type: unsigned int\n
     * Applies to: compiler and linker
     */
    CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES = 6,

    /**
     * Level of optimizations to apply to generated code (0 - 4), with 4
     * being the default and highest level of optimizations.\n
     * Option type: unsigned int\n
     * Applies to: compiler only
     */
    CU_JIT_OPTIMIZATION_LEVEL = 7,

    /**
     * No option value required. Determines the target based on the current
     * attached context (default)\n
     * Option type: No option value needed\n
     * Applies to: compiler and linker
     */
    CU_JIT_TARGET_FROM_CUCONTEXT = 8,

    /**
     * Target is chosen based on supplied ::CUjit_target.  Cannot be
     * combined with ::CU_JIT_THREADS_PER_BLOCK.\n
     * Option type: unsigned int for enumerated type ::CUjit_target\n
     * Applies to: compiler and linker
     */
    CU_JIT_TARGET = 9,

    /**
     * Specifies choice of fallback strategy if matching cubin is not found.
     * Choice is based on supplied ::CUjit_fallback.  This option cannot be
     * used with cuLink* APIs as the linker requires exact matches.\n
     * Option type: unsigned int for enumerated type ::CUjit_fallback\n
     * Applies to: compiler only
     */
    CU_JIT_FALLBACK_STRATEGY = 10,

    /**
     * Specifies whether to create debug information in output (-g)
     * (0: false, default)\n
     * Option type: int\n
     * Applies to: compiler and linker
     */
    CU_JIT_GENERATE_DEBUG_INFO = 11,

    /**
     * Generate verbose log messages (0: false, default)\n
     * Option type: int\n
     * Applies to: compiler and linker
     */
    CU_JIT_LOG_VERBOSE = 12,

    /**
     * Generate line number information (-lineinfo) (0: false, default)\n
     * Option type: int\n
     * Applies to: compiler only
     */
    CU_JIT_GENERATE_LINE_INFO = 13,

    /**
     * Specifies whether to enable caching explicitly (-dlcm) \n
     * Choice is based on supplied ::CUjit_cacheMode_enum.\n
     * Option type: unsigned int for enumerated type ::CUjit_cacheMode_enum\n
     * Applies to: compiler only
     */
    CU_JIT_CACHE_MODE = 14,

    CU_JIT_NUM_OPTIONS = 15
}

alias CUjit_option = CUjit_option_enum;

/**
 * Online compilation targets
 */
enum CUjit_target_enum
{
    CU_TARGET_COMPUTE_10 = 10, /**< Compute device class 1.0 */
    CU_TARGET_COMPUTE_11 = 11, /**< Compute device class 1.1 */
    CU_TARGET_COMPUTE_12 = 12, /**< Compute device class 1.2 */
    CU_TARGET_COMPUTE_13 = 13, /**< Compute device class 1.3 */
    CU_TARGET_COMPUTE_20 = 20, /**< Compute device class 2.0 */
    CU_TARGET_COMPUTE_21 = 21, /**< Compute device class 2.1 */
    CU_TARGET_COMPUTE_30 = 30, /**< Compute device class 3.0 */
    CU_TARGET_COMPUTE_32 = 32, /**< Compute device class 3.2 */
    CU_TARGET_COMPUTE_35 = 35, /**< Compute device class 3.5 */
    CU_TARGET_COMPUTE_37 = 37, /**< Compute device class 3.7 */
    CU_TARGET_COMPUTE_50 = 50, /**< Compute device class 5.0 */
    CU_TARGET_COMPUTE_52 = 52 /**< Compute device class 5.2 */
}

alias CUjit_target = CUjit_target_enum;

/**
 * Cubin matching fallback strategies
 */
enum CUjit_fallback_enum
{
    CU_PREFER_PTX = 0, /**< Prefer to compile ptx if exact binary match not found */

    CU_PREFER_BINARY = 1 /**< Prefer to fall back to compatible binary code if exact match not found */
}

alias CUjit_fallback = CUjit_fallback_enum;

/**
 * Caching modes for dlcm
 */
enum CUjit_cacheMode_enum
{
    CU_JIT_CACHE_OPTION_NONE = 0, /**< Compile with no -dlcm flag specified */
    CU_JIT_CACHE_OPTION_CG = 1, /**< Compile with L1 cache disabled */
    CU_JIT_CACHE_OPTION_CA = 2 /**< Compile with L1 cache enabled */
}

alias CUjit_cacheMode = CUjit_cacheMode_enum;

/**
 * Device code formats
 */
enum CUjitInputType_enum
{
    /**
     * Compiled device-class-specific device code\n
     * Applicable options: none
     */
    CU_JIT_INPUT_CUBIN = 0,

    /**
     * PTX source code\n
     * Applicable options: PTX compiler options
     */
    CU_JIT_INPUT_PTX = 1,

    /**
     * Bundle of multiple cubins and/or PTX of some device code\n
     * Applicable options: PTX compiler options, ::CU_JIT_FALLBACK_STRATEGY
     */
    CU_JIT_INPUT_FATBINARY = 2,

    /**
     * Host object with embedded device code\n
     * Applicable options: PTX compiler options, ::CU_JIT_FALLBACK_STRATEGY
     */
    CU_JIT_INPUT_OBJECT = 3,

    /**
     * Archive of host objects with embedded device code\n
     * Applicable options: PTX compiler options, ::CU_JIT_FALLBACK_STRATEGY
     */
    CU_JIT_INPUT_LIBRARY = 4,

    CU_JIT_NUM_INPUT_TYPES = 5
}

alias CUjitInputType = CUjitInputType_enum;

struct CUlinkState_st;
alias CUlinkState = CUlinkState_st*;
/* __CUDA_API_VERSION >= 5050 */

/**
 * Flags to register a graphics resource
 */
enum CUgraphicsRegisterFlags_enum
{
    CU_GRAPHICS_REGISTER_FLAGS_NONE = 0,
    CU_GRAPHICS_REGISTER_FLAGS_READ_ONLY = 1,
    CU_GRAPHICS_REGISTER_FLAGS_WRITE_DISCARD = 2,
    CU_GRAPHICS_REGISTER_FLAGS_SURFACE_LDST = 4,
    CU_GRAPHICS_REGISTER_FLAGS_TEXTURE_GATHER = 8
}

alias CUgraphicsRegisterFlags = CUgraphicsRegisterFlags_enum;

/**
 * Flags for mapping and unmapping interop resources
 */
enum CUgraphicsMapResourceFlags_enum
{
    CU_GRAPHICS_MAP_RESOURCE_FLAGS_NONE = 0,
    CU_GRAPHICS_MAP_RESOURCE_FLAGS_READ_ONLY = 1,
    CU_GRAPHICS_MAP_RESOURCE_FLAGS_WRITE_DISCARD = 2
}

alias CUgraphicsMapResourceFlags = CUgraphicsMapResourceFlags_enum;

/**
 * Array indices for cube faces
 */
enum CUarray_cubemap_face_enum
{
    CU_CUBEMAP_FACE_POSITIVE_X = 0, /**< Positive X face of cubemap */
    CU_CUBEMAP_FACE_NEGATIVE_X = 1, /**< Negative X face of cubemap */
    CU_CUBEMAP_FACE_POSITIVE_Y = 2, /**< Positive Y face of cubemap */
    CU_CUBEMAP_FACE_NEGATIVE_Y = 3, /**< Negative Y face of cubemap */
    CU_CUBEMAP_FACE_POSITIVE_Z = 4, /**< Positive Z face of cubemap */
    CU_CUBEMAP_FACE_NEGATIVE_Z = 5 /**< Negative Z face of cubemap */
}

alias CUarray_cubemap_face = CUarray_cubemap_face_enum;

/**
 * Limits
 */
enum CUlimit_enum
{
    CU_LIMIT_STACK_SIZE = 0, /**< GPU thread stack size */
    CU_LIMIT_PRINTF_FIFO_SIZE = 1, /**< GPU printf FIFO size */
    CU_LIMIT_MALLOC_HEAP_SIZE = 2, /**< GPU malloc heap size */
    CU_LIMIT_DEV_RUNTIME_SYNC_DEPTH = 3, /**< GPU device runtime launch synchronize depth */
    CU_LIMIT_DEV_RUNTIME_PENDING_LAUNCH_COUNT = 4, /**< GPU device runtime pending launch count */
    CU_LIMIT_MAX = 5
}

alias CUlimit = CUlimit_enum;

/**
 * Resource types
 */
enum CUresourcetype_enum
{
    CU_RESOURCE_TYPE_ARRAY = 0, /**< Array resoure */
    CU_RESOURCE_TYPE_MIPMAPPED_ARRAY = 1, /**< Mipmapped array resource */
    CU_RESOURCE_TYPE_LINEAR = 2, /**< Linear resource */
    CU_RESOURCE_TYPE_PITCH2D = 3 /**< Pitch 2D resource */
}

alias CUresourcetype = CUresourcetype_enum;

/**
 * Error codes
 */
enum cudaError_enum
{
    /**
     * The API call returned with no errors. In the case of query calls, this
     * can also mean that the operation being queried is complete (see
     * ::cuEventQuery() and ::cuStreamQuery()).
     */
    CUDA_SUCCESS = 0,

    /**
     * This indicates that one or more of the parameters passed to the API call
     * is not within an acceptable range of values.
     */
    CUDA_ERROR_INVALID_VALUE = 1,

    /**
     * The API call failed because it was unable to allocate enough memory to
     * perform the requested operation.
     */
    CUDA_ERROR_OUT_OF_MEMORY = 2,

    /**
     * This indicates that the CUDA driver has not been initialized with
     * ::cuInit() or that initialization has failed.
     */
    CUDA_ERROR_NOT_INITIALIZED = 3,

    /**
     * This indicates that the CUDA driver is in the process of shutting down.
     */
    CUDA_ERROR_DEINITIALIZED = 4,

    /**
     * This indicates profiler is not initialized for this run. This can
     * happen when the application is running with external profiling tools
     * like visual profiler.
     */
    CUDA_ERROR_PROFILER_DISABLED = 5,

    /**
     * \deprecated
     * This error return is deprecated as of CUDA 5.0. It is no longer an error
     * to attempt to enable/disable the profiling via ::cuProfilerStart or
     * ::cuProfilerStop without initialization.
     */
    CUDA_ERROR_PROFILER_NOT_INITIALIZED = 6,

    /**
     * \deprecated
     * This error return is deprecated as of CUDA 5.0. It is no longer an error
     * to call cuProfilerStart() when profiling is already enabled.
     */
    CUDA_ERROR_PROFILER_ALREADY_STARTED = 7,

    /**
     * \deprecated
     * This error return is deprecated as of CUDA 5.0. It is no longer an error
     * to call cuProfilerStop() when profiling is already disabled.
     */
    CUDA_ERROR_PROFILER_ALREADY_STOPPED = 8,

    /**
     * This indicates that no CUDA-capable devices were detected by the installed
     * CUDA driver.
     */
    CUDA_ERROR_NO_DEVICE = 100,

    /**
     * This indicates that the device ordinal supplied by the user does not
     * correspond to a valid CUDA device.
     */
    CUDA_ERROR_INVALID_DEVICE = 101,

    /**
     * This indicates that the device kernel image is invalid. This can also
     * indicate an invalid CUDA module.
     */
    CUDA_ERROR_INVALID_IMAGE = 200,

    /**
     * This most frequently indicates that there is no context bound to the
     * current thread. This can also be returned if the context passed to an
     * API call is not a valid handle (such as a context that has had
     * ::cuCtxDestroy() invoked on it). This can also be returned if a user
     * mixes different API versions (i.e. 3010 context with 3020 API calls).
     * See ::cuCtxGetApiVersion() for more details.
     */
    CUDA_ERROR_INVALID_CONTEXT = 201,

    /**
     * This indicated that the context being supplied as a parameter to the
     * API call was already the active context.
     * \deprecated
     * This error return is deprecated as of CUDA 3.2. It is no longer an
     * error to attempt to push the active context via ::cuCtxPushCurrent().
     */
    CUDA_ERROR_CONTEXT_ALREADY_CURRENT = 202,

    /**
     * This indicates that a map or register operation has failed.
     */
    CUDA_ERROR_MAP_FAILED = 205,

    /**
     * This indicates that an unmap or unregister operation has failed.
     */
    CUDA_ERROR_UNMAP_FAILED = 206,

    /**
     * This indicates that the specified array is currently mapped and thus
     * cannot be destroyed.
     */
    CUDA_ERROR_ARRAY_IS_MAPPED = 207,

    /**
     * This indicates that the resource is already mapped.
     */
    CUDA_ERROR_ALREADY_MAPPED = 208,

    /**
     * This indicates that there is no kernel image available that is suitable
     * for the device. This can occur when a user specifies code generation
     * options for a particular CUDA source file that do not include the
     * corresponding device configuration.
     */
    CUDA_ERROR_NO_BINARY_FOR_GPU = 209,

    /**
     * This indicates that a resource has already been acquired.
     */
    CUDA_ERROR_ALREADY_ACQUIRED = 210,

    /**
     * This indicates that a resource is not mapped.
     */
    CUDA_ERROR_NOT_MAPPED = 211,

    /**
     * This indicates that a mapped resource is not available for access as an
     * array.
     */
    CUDA_ERROR_NOT_MAPPED_AS_ARRAY = 212,

    /**
     * This indicates that a mapped resource is not available for access as a
     * pointer.
     */
    CUDA_ERROR_NOT_MAPPED_AS_POINTER = 213,

    /**
     * This indicates that an uncorrectable ECC error was detected during
     * execution.
     */
    CUDA_ERROR_ECC_UNCORRECTABLE = 214,

    /**
     * This indicates that the ::CUlimit passed to the API call is not
     * supported by the active device.
     */
    CUDA_ERROR_UNSUPPORTED_LIMIT = 215,

    /**
     * This indicates that the ::CUcontext passed to the API call can
     * only be bound to a single CPU thread at a time but is already
     * bound to a CPU thread.
     */
    CUDA_ERROR_CONTEXT_ALREADY_IN_USE = 216,

    /**
     * This indicates that peer access is not supported across the given
     * devices.
     */
    CUDA_ERROR_PEER_ACCESS_UNSUPPORTED = 217,

    /**
     * This indicates that a PTX JIT compilation failed.
     */
    CUDA_ERROR_INVALID_PTX = 218,

    /**
     * This indicates an error with OpenGL or DirectX context.
     */
    CUDA_ERROR_INVALID_GRAPHICS_CONTEXT = 219,

    /**
     * This indicates that the device kernel source is invalid.
     */
    CUDA_ERROR_INVALID_SOURCE = 300,

    /**
     * This indicates that the file specified was not found.
     */
    CUDA_ERROR_FILE_NOT_FOUND = 301,

    /**
     * This indicates that a link to a shared object failed to resolve.
     */
    CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND = 302,

    /**
     * This indicates that initialization of a shared object failed.
     */
    CUDA_ERROR_SHARED_OBJECT_INIT_FAILED = 303,

    /**
     * This indicates that an OS call failed.
     */
    CUDA_ERROR_OPERATING_SYSTEM = 304,

    /**
     * This indicates that a resource handle passed to the API call was not
     * valid. Resource handles are opaque types like ::CUstream and ::CUevent.
     */
    CUDA_ERROR_INVALID_HANDLE = 400,

    /**
     * This indicates that a named symbol was not found. Examples of symbols
     * are global/constant variable names, texture names, and surface names.
     */
    CUDA_ERROR_NOT_FOUND = 500,

    /**
     * This indicates that asynchronous operations issued previously have not
     * completed yet. This result is not actually an error, but must be indicated
     * differently than ::CUDA_SUCCESS (which indicates completion). Calls that
     * may return this value include ::cuEventQuery() and ::cuStreamQuery().
     */
    CUDA_ERROR_NOT_READY = 600,

    /**
     * While executing a kernel, the device encountered a
     * load or store instruction on an invalid memory address.
     * The context cannot be used, so it must be destroyed (and a new one should be created).
     * All existing device memory allocations from this context are invalid
     * and must be reconstructed if the program is to continue using CUDA.
     */
    CUDA_ERROR_ILLEGAL_ADDRESS = 700,

    /**
     * This indicates that a launch did not occur because it did not have
     * appropriate resources. This error usually indicates that the user has
     * attempted to pass too many arguments to the device kernel, or the
     * kernel launch specifies too many threads for the kernel's register
     * count. Passing arguments of the wrong size (i.e. a 64-bit pointer
     * when a 32-bit int is expected) is equivalent to passing too many
     * arguments and can also result in this error.
     */
    CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES = 701,

    /**
     * This indicates that the device kernel took too long to execute. This can
     * only occur if timeouts are enabled - see the device attribute
     * ::CU_DEVICE_ATTRIBUTE_KERNEL_EXEC_TIMEOUT for more information. The
     * context cannot be used (and must be destroyed similar to
     * ::CUDA_ERROR_LAUNCH_FAILED). All existing device memory allocations from
     * this context are invalid and must be reconstructed if the program is to
     * continue using CUDA.
     */
    CUDA_ERROR_LAUNCH_TIMEOUT = 702,

    /**
     * This error indicates a kernel launch that uses an incompatible texturing
     * mode.
     */
    CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING = 703,

    /**
     * This error indicates that a call to ::cuCtxEnablePeerAccess() is
     * trying to re-enable peer access to a context which has already
     * had peer access to it enabled.
     */
    CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED = 704,

    /**
     * This error indicates that ::cuCtxDisablePeerAccess() is
     * trying to disable peer access which has not been enabled yet
     * via ::cuCtxEnablePeerAccess().
     */
    CUDA_ERROR_PEER_ACCESS_NOT_ENABLED = 705,

    /**
     * This error indicates that the primary context for the specified device
     * has already been initialized.
     */
    CUDA_ERROR_PRIMARY_CONTEXT_ACTIVE = 708,

    /**
     * This error indicates that the context current to the calling thread
     * has been destroyed using ::cuCtxDestroy, or is a primary context which
     * has not yet been initialized.
     */
    CUDA_ERROR_CONTEXT_IS_DESTROYED = 709,

    /**
     * A device-side assert triggered during kernel execution. The context
     * cannot be used anymore, and must be destroyed. All existing device
     * memory allocations from this context are invalid and must be
     * reconstructed if the program is to continue using CUDA.
     */
    CUDA_ERROR_ASSERT = 710,

    /**
     * This error indicates that the hardware resources required to enable
     * peer access have been exhausted for one or more of the devices
     * passed to ::cuCtxEnablePeerAccess().
     */
    CUDA_ERROR_TOO_MANY_PEERS = 711,

    /**
     * This error indicates that the memory range passed to ::cuMemHostRegister()
     * has already been registered.
     */
    CUDA_ERROR_HOST_MEMORY_ALREADY_REGISTERED = 712,

    /**
     * This error indicates that the pointer passed to ::cuMemHostUnregister()
     * does not correspond to any currently registered memory region.
     */
    CUDA_ERROR_HOST_MEMORY_NOT_REGISTERED = 713,

    /**
     * While executing a kernel, the device encountered a stack error.
     * This can be due to stack corruption or exceeding the stack size limit.
     * The context cannot be used, so it must be destroyed (and a new one should be created).
     * All existing device memory allocations from this context are invalid
     * and must be reconstructed if the program is to continue using CUDA.
     */
    CUDA_ERROR_HARDWARE_STACK_ERROR = 714,

    /**
     * While executing a kernel, the device encountered an illegal instruction.
     * The context cannot be used, so it must be destroyed (and a new one should be created).
     * All existing device memory allocations from this context are invalid
     * and must be reconstructed if the program is to continue using CUDA.
     */
    CUDA_ERROR_ILLEGAL_INSTRUCTION = 715,

    /**
     * While executing a kernel, the device encountered a load or store instruction
     * on a memory address which is not aligned.
     * The context cannot be used, so it must be destroyed (and a new one should be created).
     * All existing device memory allocations from this context are invalid
     * and must be reconstructed if the program is to continue using CUDA.
     */
    CUDA_ERROR_MISALIGNED_ADDRESS = 716,

    /**
     * While executing a kernel, the device encountered an instruction
     * which can only operate on memory locations in certain address spaces
     * (global, shared, or local), but was supplied a memory address not
     * belonging to an allowed address space.
     * The context cannot be used, so it must be destroyed (and a new one should be created).
     * All existing device memory allocations from this context are invalid
     * and must be reconstructed if the program is to continue using CUDA.
     */
    CUDA_ERROR_INVALID_ADDRESS_SPACE = 717,

    /**
     * While executing a kernel, the device program counter wrapped its address space.
     * The context cannot be used, so it must be destroyed (and a new one should be created).
     * All existing device memory allocations from this context are invalid
     * and must be reconstructed if the program is to continue using CUDA.
     */
    CUDA_ERROR_INVALID_PC = 718,

    /**
     * An exception occurred on the device while executing a kernel. Common
     * causes include dereferencing an invalid device pointer and accessing
     * out of bounds shared memory. The context cannot be used, so it must
     * be destroyed (and a new one should be created). All existing device
     * memory allocations from this context are invalid and must be
     * reconstructed if the program is to continue using CUDA.
     */
    CUDA_ERROR_LAUNCH_FAILED = 719,

    /**
     * This error indicates that the attempted operation is not permitted.
     */
    CUDA_ERROR_NOT_PERMITTED = 800,

    /**
     * This error indicates that the attempted operation is not supported
     * on the current system or device.
     */
    CUDA_ERROR_NOT_SUPPORTED = 801,

    /**
     * This indicates that an unknown internal error has occurred.
     */
    CUDA_ERROR_UNKNOWN = 999
}

alias CUresult = cudaError_enum;

alias CUstreamCallback = void function (CUstream hStream, CUresult status, void* userData);

alias CUoccupancyB2DSize = c_ulong function (int blockSize);

enum CU_MEMHOSTALLOC_PORTABLE = 0x01;

enum CU_MEMHOSTALLOC_DEVICEMAP = 0x02;

enum CU_MEMHOSTALLOC_WRITECOMBINED = 0x04;

enum CU_MEMHOSTREGISTER_PORTABLE = 0x01;

enum CU_MEMHOSTREGISTER_DEVICEMAP = 0x02;

enum CU_MEMHOSTREGISTER_IOMEMORY = 0x04;

struct CUDA_MEMCPY2D_st
{
    size_t srcXInBytes; /**< Source X in bytes */
    size_t srcY; /**< Source Y */

    CUmemorytype srcMemoryType; /**< Source memory type (host, device, array) */
    const(void)* srcHost; /**< Source host pointer */
    CUdeviceptr srcDevice; /**< Source device pointer */
    CUarray srcArray; /**< Source array reference */
    size_t srcPitch; /**< Source pitch (ignored when src is array) */

    size_t dstXInBytes; /**< Destination X in bytes */
    size_t dstY; /**< Destination Y */

    CUmemorytype dstMemoryType; /**< Destination memory type (host, device, array) */
    void* dstHost; /**< Destination host pointer */
    CUdeviceptr dstDevice; /**< Destination device pointer */
    CUarray dstArray; /**< Destination array reference */
    size_t dstPitch; /**< Destination pitch (ignored when dst is array) */

    size_t WidthInBytes; /**< Width of 2D memory copy in bytes */
    size_t Height; /**< Height of 2D memory copy */
}

alias CUDA_MEMCPY2D = CUDA_MEMCPY2D_st;

/**
 * 3D memory copy parameters
 */
struct CUDA_MEMCPY3D_st
{
    size_t srcXInBytes; /**< Source X in bytes */
    size_t srcY; /**< Source Y */
    size_t srcZ; /**< Source Z */
    size_t srcLOD; /**< Source LOD */
    CUmemorytype srcMemoryType; /**< Source memory type (host, device, array) */
    const(void)* srcHost; /**< Source host pointer */
    CUdeviceptr srcDevice; /**< Source device pointer */
    CUarray srcArray; /**< Source array reference */
    void* reserved0; /**< Must be NULL */
    size_t srcPitch; /**< Source pitch (ignored when src is array) */
    size_t srcHeight; /**< Source height (ignored when src is array; may be 0 if Depth==1) */

    size_t dstXInBytes; /**< Destination X in bytes */
    size_t dstY; /**< Destination Y */
    size_t dstZ; /**< Destination Z */
    size_t dstLOD; /**< Destination LOD */
    CUmemorytype dstMemoryType; /**< Destination memory type (host, device, array) */
    void* dstHost; /**< Destination host pointer */
    CUdeviceptr dstDevice; /**< Destination device pointer */
    CUarray dstArray; /**< Destination array reference */
    void* reserved1; /**< Must be NULL */
    size_t dstPitch; /**< Destination pitch (ignored when dst is array) */
    size_t dstHeight; /**< Destination height (ignored when dst is array; may be 0 if Depth==1) */

    size_t WidthInBytes; /**< Width of 3D memory copy in bytes */
    size_t Height; /**< Height of 3D memory copy */
    size_t Depth; /**< Depth of 3D memory copy */
}

alias CUDA_MEMCPY3D = CUDA_MEMCPY3D_st;

/**
 * 3D memory cross-context copy parameters
 */
struct CUDA_MEMCPY3D_PEER_st
{
    size_t srcXInBytes; /**< Source X in bytes */
    size_t srcY; /**< Source Y */
    size_t srcZ; /**< Source Z */
    size_t srcLOD; /**< Source LOD */
    CUmemorytype srcMemoryType; /**< Source memory type (host, device, array) */
    const(void)* srcHost; /**< Source host pointer */
    CUdeviceptr srcDevice; /**< Source device pointer */
    CUarray srcArray; /**< Source array reference */
    CUcontext srcContext; /**< Source context (ignored with srcMemoryType is ::CU_MEMORYTYPE_ARRAY) */
    size_t srcPitch; /**< Source pitch (ignored when src is array) */
    size_t srcHeight; /**< Source height (ignored when src is array; may be 0 if Depth==1) */

    size_t dstXInBytes; /**< Destination X in bytes */
    size_t dstY; /**< Destination Y */
    size_t dstZ; /**< Destination Z */
    size_t dstLOD; /**< Destination LOD */
    CUmemorytype dstMemoryType; /**< Destination memory type (host, device, array) */
    void* dstHost; /**< Destination host pointer */
    CUdeviceptr dstDevice; /**< Destination device pointer */
    CUarray dstArray; /**< Destination array reference */
    CUcontext dstContext; /**< Destination context (ignored with dstMemoryType is ::CU_MEMORYTYPE_ARRAY) */
    size_t dstPitch; /**< Destination pitch (ignored when dst is array) */
    size_t dstHeight; /**< Destination height (ignored when dst is array; may be 0 if Depth==1) */

    size_t WidthInBytes; /**< Width of 3D memory copy in bytes */
    size_t Height; /**< Height of 3D memory copy */
    size_t Depth; /**< Depth of 3D memory copy */
}

alias CUDA_MEMCPY3D_PEER = CUDA_MEMCPY3D_PEER_st;

/**
 * Array descriptor
 */
struct CUDA_ARRAY_DESCRIPTOR_st
{
    size_t Width; /**< Width of array */
    size_t Height; /**< Height of array */

    CUarray_format Format; /**< Array format */
    uint NumChannels; /**< Channels per array element */
}

alias CUDA_ARRAY_DESCRIPTOR = CUDA_ARRAY_DESCRIPTOR_st;

/**
 * 3D array descriptor
 */
struct CUDA_ARRAY3D_DESCRIPTOR_st
{
    size_t Width; /**< Width of 3D array */
    size_t Height; /**< Height of 3D array */
    size_t Depth; /**< Depth of 3D array */

    CUarray_format Format; /**< Array format */
    uint NumChannels; /**< Channels per array element */
    uint Flags; /**< Flags */
}

alias CUDA_ARRAY3D_DESCRIPTOR = CUDA_ARRAY3D_DESCRIPTOR_st;

/* __CUDA_API_VERSION >= 3020 */

/**
 * CUDA Resource descriptor
 */
struct CUDA_RESOURCE_DESC_st
{
    CUresourcetype resType; /**< Resource type */

    /**< CUDA array */

    /**< CUDA mipmapped array */

    /**< Device pointer */
    /**< Array format */
    /**< Channels per array element */
    /**< Size in bytes */

    /**< Device pointer */
    /**< Array format */
    /**< Channels per array element */
    /**< Width of the array in elements */
    /**< Height of the array in elements */
    /**< Pitch between two rows in bytes */
    union _Anonymous_0
    {
        struct _Anonymous_1
        {
            CUarray hArray;
        }

        _Anonymous_1 array;

        struct _Anonymous_2
        {
            CUmipmappedArray hMipmappedArray;
        }

        _Anonymous_2 mipmap;

        struct _Anonymous_3
        {
            CUdeviceptr devPtr;
            CUarray_format format;
            uint numChannels;
            size_t sizeInBytes;
        }

        _Anonymous_3 linear;

        struct _Anonymous_4
        {
            CUdeviceptr devPtr;
            CUarray_format format;
            uint numChannels;
            size_t width;
            size_t height;
            size_t pitchInBytes;
        }

        _Anonymous_4 pitch2D;

        struct _Anonymous_5
        {
            int[32] reserved;
        }

        _Anonymous_5 reserved;
    }

    _Anonymous_0 res;

    uint flags; /**< Flags (must be zero) */
}

alias CUDA_RESOURCE_DESC = CUDA_RESOURCE_DESC_st;

/**
 * Texture descriptor
 */
struct CUDA_TEXTURE_DESC_st
{
    CUaddress_mode[3] addressMode; /**< Address modes */
    CUfilter_mode filterMode; /**< Filter mode */
    uint flags; /**< Flags */
    uint maxAnisotropy; /**< Maximum anisotropy ratio */
    CUfilter_mode mipmapFilterMode; /**< Mipmap filter mode */
    float mipmapLevelBias; /**< Mipmap level bias */
    float minMipmapLevelClamp; /**< Mipmap minimum level clamp */
    float maxMipmapLevelClamp; /**< Mipmap maximum level clamp */
    int[16] reserved;
}

alias CUDA_TEXTURE_DESC = CUDA_TEXTURE_DESC_st;

/**
 * Resource view format
 */
enum CUresourceViewFormat_enum
{
    CU_RES_VIEW_FORMAT_NONE = 0, /**< No resource view format (use underlying resource format) */
    CU_RES_VIEW_FORMAT_UINT_1X8 = 1, /**< 1 channel unsigned 8-bit integers */
    CU_RES_VIEW_FORMAT_UINT_2X8 = 2, /**< 2 channel unsigned 8-bit integers */
    CU_RES_VIEW_FORMAT_UINT_4X8 = 3, /**< 4 channel unsigned 8-bit integers */
    CU_RES_VIEW_FORMAT_SINT_1X8 = 4, /**< 1 channel signed 8-bit integers */
    CU_RES_VIEW_FORMAT_SINT_2X8 = 5, /**< 2 channel signed 8-bit integers */
    CU_RES_VIEW_FORMAT_SINT_4X8 = 6, /**< 4 channel signed 8-bit integers */
    CU_RES_VIEW_FORMAT_UINT_1X16 = 7, /**< 1 channel unsigned 16-bit integers */
    CU_RES_VIEW_FORMAT_UINT_2X16 = 8, /**< 2 channel unsigned 16-bit integers */
    CU_RES_VIEW_FORMAT_UINT_4X16 = 9, /**< 4 channel unsigned 16-bit integers */
    CU_RES_VIEW_FORMAT_SINT_1X16 = 10, /**< 1 channel signed 16-bit integers */
    CU_RES_VIEW_FORMAT_SINT_2X16 = 11, /**< 2 channel signed 16-bit integers */
    CU_RES_VIEW_FORMAT_SINT_4X16 = 12, /**< 4 channel signed 16-bit integers */
    CU_RES_VIEW_FORMAT_UINT_1X32 = 13, /**< 1 channel unsigned 32-bit integers */
    CU_RES_VIEW_FORMAT_UINT_2X32 = 14, /**< 2 channel unsigned 32-bit integers */
    CU_RES_VIEW_FORMAT_UINT_4X32 = 15, /**< 4 channel unsigned 32-bit integers */
    CU_RES_VIEW_FORMAT_SINT_1X32 = 16, /**< 1 channel signed 32-bit integers */
    CU_RES_VIEW_FORMAT_SINT_2X32 = 17, /**< 2 channel signed 32-bit integers */
    CU_RES_VIEW_FORMAT_SINT_4X32 = 18, /**< 4 channel signed 32-bit integers */
    CU_RES_VIEW_FORMAT_FLOAT_1X16 = 19, /**< 1 channel 16-bit floating point */
    CU_RES_VIEW_FORMAT_FLOAT_2X16 = 20, /**< 2 channel 16-bit floating point */
    CU_RES_VIEW_FORMAT_FLOAT_4X16 = 21, /**< 4 channel 16-bit floating point */
    CU_RES_VIEW_FORMAT_FLOAT_1X32 = 22, /**< 1 channel 32-bit floating point */
    CU_RES_VIEW_FORMAT_FLOAT_2X32 = 23, /**< 2 channel 32-bit floating point */
    CU_RES_VIEW_FORMAT_FLOAT_4X32 = 24, /**< 4 channel 32-bit floating point */
    CU_RES_VIEW_FORMAT_UNSIGNED_BC1 = 25, /**< Block compressed 1 */
    CU_RES_VIEW_FORMAT_UNSIGNED_BC2 = 26, /**< Block compressed 2 */
    CU_RES_VIEW_FORMAT_UNSIGNED_BC3 = 27, /**< Block compressed 3 */
    CU_RES_VIEW_FORMAT_UNSIGNED_BC4 = 28, /**< Block compressed 4 unsigned */
    CU_RES_VIEW_FORMAT_SIGNED_BC4 = 29, /**< Block compressed 4 signed */
    CU_RES_VIEW_FORMAT_UNSIGNED_BC5 = 30, /**< Block compressed 5 unsigned */
    CU_RES_VIEW_FORMAT_SIGNED_BC5 = 31, /**< Block compressed 5 signed */
    CU_RES_VIEW_FORMAT_UNSIGNED_BC6H = 32, /**< Block compressed 6 unsigned half-float */
    CU_RES_VIEW_FORMAT_SIGNED_BC6H = 33, /**< Block compressed 6 signed half-float */
    CU_RES_VIEW_FORMAT_UNSIGNED_BC7 = 34 /**< Block compressed 7 */
}

alias CUresourceViewFormat = CUresourceViewFormat_enum;

/**
 * Resource view descriptor
 */
struct CUDA_RESOURCE_VIEW_DESC_st
{
    CUresourceViewFormat format; /**< Resource view format */
    size_t width; /**< Width of the resource view */
    size_t height; /**< Height of the resource view */
    size_t depth; /**< Depth of the resource view */
    uint firstMipmapLevel; /**< First defined mipmap level */
    uint lastMipmapLevel; /**< Last defined mipmap level */
    uint firstLayer; /**< First layer index */
    uint lastLayer; /**< Last layer index */
    uint[16] reserved;
}

alias CUDA_RESOURCE_VIEW_DESC = CUDA_RESOURCE_VIEW_DESC_st;

/**
 * GPU Direct v3 tokens
 */
struct CUDA_POINTER_ATTRIBUTE_P2P_TOKENS_st
{
    ulong p2pToken;
    uint vaSpaceToken;
}

alias CUDA_POINTER_ATTRIBUTE_P2P_TOKENS = CUDA_POINTER_ATTRIBUTE_P2P_TOKENS_st;

/* __CUDA_API_VERSION >= 5000 */

/**
 * If set, the CUDA array is a collection of layers, where each layer is either a 1D
 * or a 2D array and the Depth member of CUDA_ARRAY3D_DESCRIPTOR specifies the number
 * of layers, not the depth of a 3D array.
 */
enum CUDA_ARRAY3D_LAYERED = 0x01;

/**
 * Deprecated, use CUDA_ARRAY3D_LAYERED
 */
enum CUDA_ARRAY3D_2DARRAY = 0x01;

/**
 * This flag must be set in order to bind a surface reference
 * to the CUDA array
 */
enum CUDA_ARRAY3D_SURFACE_LDST = 0x02;

/**
 * If set, the CUDA array is a collection of six 2D arrays, representing faces of a cube. The
 * width of such a CUDA array must be equal to its height, and Depth must be six.
 * If ::CUDA_ARRAY3D_LAYERED flag is also set, then the CUDA array is a collection of cubemaps
 * and Depth must be a multiple of six.
 */
enum CUDA_ARRAY3D_CUBEMAP = 0x04;

/**
 * This flag must be set in order to perform texture gather operations
 * on a CUDA array.
 */
enum CUDA_ARRAY3D_TEXTURE_GATHER = 0x08;

/**
 * This flag if set indicates that the CUDA
 * array is a DEPTH_TEXTURE.
*/
enum CUDA_ARRAY3D_DEPTH_TEXTURE = 0x10;

/**
 * Override the texref format with a format inferred from the array.
 * Flag for ::cuTexRefSetArray()
 */
enum CU_TRSA_OVERRIDE_FORMAT = 0x01;

/**
 * Read the texture as integers rather than promoting the values to floats
 * in the range [0,1].
 * Flag for ::cuTexRefSetFlags()
 */
enum CU_TRSF_READ_AS_INTEGER = 0x01;

/**
 * Use normalized texture coordinates in the range [0,1) instead of [0,dim).
 * Flag for ::cuTexRefSetFlags()
 */
enum CU_TRSF_NORMALIZED_COORDINATES = 0x02;

/**
 * Perform sRGB->linear conversion during texture read.
 * Flag for ::cuTexRefSetFlags()
 */
enum CU_TRSF_SRGB = 0x10;

/**
 * End of array terminator for the \p extra parameter to
 * ::cuLaunchKernel
 */
enum CU_LAUNCH_PARAM_END = cast(void*) 0x00;

/**
 * Indicator that the next value in the \p extra parameter to
 * ::cuLaunchKernel will be a pointer to a buffer containing all kernel
 * parameters used for launching kernel \p f.  This buffer needs to
 * honor all alignment/padding requirements of the individual parameters.
 * If ::CU_LAUNCH_PARAM_BUFFER_SIZE is not also specified in the
 * \p extra array, then ::CU_LAUNCH_PARAM_BUFFER_POINTER will have no
 * effect.
 */
enum CU_LAUNCH_PARAM_BUFFER_POINTER = cast(void*) 0x01;

/**
 * Indicator that the next value in the \p extra parameter to
 * ::cuLaunchKernel will be a pointer to a size_t which contains the
 * size of the buffer specified with ::CU_LAUNCH_PARAM_BUFFER_POINTER.
 * It is required that ::CU_LAUNCH_PARAM_BUFFER_POINTER also be specified
 * in the \p extra array if the value associated with
 * ::CU_LAUNCH_PARAM_BUFFER_SIZE is not zero.
 */
enum CU_LAUNCH_PARAM_BUFFER_SIZE = cast(void*) 0x02;

/**
 * For texture references loaded into the module, use default texunit from
 * texture reference.
 */
enum CU_PARAM_TR_DEFAULT = -1;

CUresult cuGetErrorString (CUresult error, const(char*)* pStr);

CUresult cuGetErrorName (CUresult error, const(char*)* pStr);

CUresult cuInit (uint Flags);

CUresult cuDriverGetVersion (int* driverVersion);

CUresult cuDeviceGet (CUdevice* device, int ordinal);

CUresult cuDeviceGetCount (int* count);

CUresult cuDeviceGetName (char[64]* name, int len, CUdevice dev);

CUresult cuDeviceTotalMem_v2 (size_t* bytes, CUdevice dev);

CUresult cuDeviceGetAttribute (int* pi, CUdevice_attribute attrib, CUdevice dev);

CUresult cuDeviceGetProperties (CUdevprop* prop, CUdevice dev);

CUresult cuDeviceComputeCapability (int* major, int* minor, CUdevice dev);

CUresult cuDevicePrimaryCtxRetain (CUcontext* pctx, CUdevice dev);

CUresult cuDevicePrimaryCtxRelease (CUdevice dev);

CUresult cuDevicePrimaryCtxSetFlags (CUdevice dev, uint flags);

CUresult cuDevicePrimaryCtxGetState (CUdevice dev, uint* flags, int* active);

CUresult cuDevicePrimaryCtxReset (CUdevice dev);

CUresult cuCtxCreate_v2 (CUcontext* pctx, uint flags, CUdevice dev);

CUresult cuCtxDestroy_v2 (CUcontext ctx);

CUresult cuCtxPushCurrent_v2 (CUcontext ctx);

CUresult cuCtxPopCurrent_v2 (CUcontext* pctx);

CUresult cuCtxSetCurrent (CUcontext ctx);

CUresult cuCtxGetCurrent (CUcontext* pctx);

CUresult cuCtxGetDevice (CUdevice* device);

CUresult cuCtxGetFlags (uint* flags);


CUresult cuCtxSynchronize ();

CUresult cuCtxSetLimit (CUlimit limit, size_t value);


CUresult cuCtxGetLimit (size_t* pvalue, CUlimit limit);

CUresult cuCtxGetCacheConfig (CUfunc_cache* pconfig);

CUresult cuCtxSetCacheConfig (CUfunc_cache config);

/**
 * \brief Returns the current shared memory configuration for the current context.
 *
 * This function will return in \p pConfig the current size of shared memory banks
 * in the current context. On devices with configurable shared memory banks,
 * ::cuCtxSetSharedMemConfig can be used to change this setting, so that all
 * subsequent kernel launches will by default use the new bank size. When
 * ::cuCtxGetSharedMemConfig is called on devices without configurable shared
 * memory, it will return the fixed bank size of the hardware.
 *
 * The returned bank configurations can be either:
 * - ::CU_SHARED_MEM_CONFIG_FOUR_BYTE_BANK_SIZE:  shared memory bank width is
 *   four bytes.
 * - ::CU_SHARED_MEM_CONFIG_EIGHT_BYTE_BANK_SIZE: shared memory bank width will
 *   eight bytes.
 *
 * \param pConfig - returned shared memory configuration
 * \return
 * ::CUDA_SUCCESS,
 * ::CUDA_ERROR_DEINITIALIZED,
 * ::CUDA_ERROR_NOT_INITIALIZED,
 * ::CUDA_ERROR_INVALID_CONTEXT,
 * ::CUDA_ERROR_INVALID_VALUE
 * \notefnerr
 *
 * \sa ::cuCtxCreate,
 * ::cuCtxDestroy,
 * ::cuCtxGetApiVersion,
 * ::cuCtxGetCacheConfig,
 * ::cuCtxGetDevice,
 * ::cuCtxGetFlags,
 * ::cuCtxGetLimit,
 * ::cuCtxPopCurrent,
 * ::cuCtxPushCurrent,
 * ::cuCtxSetLimit,
 * ::cuCtxSynchronize,
 * ::cuCtxGetSharedMemConfig,
 * ::cuFuncSetCacheConfig,
 */
CUresult cuCtxGetSharedMemConfig (CUsharedconfig* pConfig);

CUresult cuCtxSetSharedMemConfig (CUsharedconfig config);

/**
 * \brief Gets the context's API version.
 *
 * Returns a version number in \p version corresponding to the capabilities of
 * the context (e.g. 3010 or 3020), which library developers can use to direct
 * callers to a specific API version. If \p ctx is NULL, returns the API version
 * used to create the currently bound context.
 *
 * Note that new API versions are only introduced when context capabilities are
 * changed that break binary compatibility, so the API version and driver version
 * may be different. For example, it is valid for the API version to be 3020 while
 * the driver version is 4020.
 *
 * \param ctx     - Context to check
 * \param version - Pointer to version
 *
 * \return
 * ::CUDA_SUCCESS,
 * ::CUDA_ERROR_DEINITIALIZED,
 * ::CUDA_ERROR_NOT_INITIALIZED,
 * ::CUDA_ERROR_INVALID_CONTEXT,
 * ::CUDA_ERROR_UNKNOWN
 * \notefnerr
 *
 * \sa ::cuCtxCreate,
 * ::cuCtxDestroy,
 * ::cuCtxGetDevice,
 * ::cuCtxGetFlags,
 * ::cuCtxGetLimit,
 * ::cuCtxPopCurrent,
 * ::cuCtxPushCurrent,
 * ::cuCtxSetCacheConfig,
 * ::cuCtxSetLimit,
 * ::cuCtxSynchronize
 */
CUresult cuCtxGetApiVersion (CUcontext ctx, uint* version_);


CUresult cuCtxGetStreamPriorityRange (int* leastPriority, int* greatestPriority);


CUresult cuCtxAttach (CUcontext* pctx, uint flags);

CUresult cuCtxDetach (CUcontext ctx);

CUresult cuModuleLoad (CUmodule* module_, const(char)* fname);

CUresult cuModuleLoadData (CUmodule* module_, const(void)* image);

CUresult cuModuleLoadDataEx (CUmodule* module_, const(void)* image, uint numOptions, CUjit_option* options, void** optionValues);

CUresult cuModuleLoadFatBinary (CUmodule* module_, const(void)* fatCubin);

CUresult cuModuleUnload (CUmodule hmod);

CUresult cuModuleGetFunction (CUfunction* hfunc, CUmodule hmod, const(char)* name);

CUresult cuModuleGetGlobal_v2 (CUdeviceptr* dptr, size_t* bytes, CUmodule hmod, const(char)* name);

CUresult cuModuleGetTexRef (CUtexref* pTexRef, CUmodule hmod, const(char)* name);

CUresult cuModuleGetSurfRef (CUsurfref* pSurfRef, CUmodule hmod, const(char)* name);

CUresult cuLinkCreate_v2 (
    uint numOptions,
    CUjit_option* options,
    void** optionValues,
    CUlinkState* stateOut);


CUresult cuLinkAddData_v2 (
    CUlinkState state,
    CUjitInputType type,
    void* data,
    size_t size,
    const(char)* name,
    uint numOptions,
    CUjit_option* options,
    void** optionValues);

CUresult cuLinkAddFile_v2 (
    CUlinkState state,
    CUjitInputType type,
    const(char)* path,
    uint numOptions,
    CUjit_option* options,
    void** optionValues);

CUresult cuLinkComplete (CUlinkState state, void** cubinOut, size_t* sizeOut);

CUresult cuLinkDestroy (CUlinkState state);

CUresult cuMemGetInfo_v2 (size_t* free, size_t* total);

CUresult cuMemAlloc_v2 (CUdeviceptr* dptr, size_t bytesize);

CUresult cuMemAllocPitch_v2 (CUdeviceptr* dptr, size_t* pPitch, size_t WidthInBytes, size_t Height, uint ElementSizeBytes);

CUresult cuMemFree_v2 (CUdeviceptr dptr);

CUresult cuMemGetAddressRange_v2 (CUdeviceptr* pbase, size_t* psize, CUdeviceptr dptr);

CUresult cuMemAllocHost_v2 (void** pp, size_t bytesize);

CUresult cuMemFreeHost (void* p);

CUresult cuMemHostAlloc (void** pp, size_t bytesize, uint Flags);

CUresult cuMemHostGetDevicePointer_v2 (CUdeviceptr* pdptr, void* p, uint Flags);

CUresult cuMemHostGetFlags (uint* pFlags, void* p);

CUresult cuMemAllocManaged (CUdeviceptr* dptr, size_t bytesize, uint flags);

CUresult cuDeviceGetByPCIBusId (CUdevice* dev, const(char)* pciBusId);

CUresult cuDeviceGetPCIBusId (char* pciBusId, int len, CUdevice dev);


CUresult cuIpcGetEventHandle (CUipcEventHandle* pHandle, CUevent event);


CUresult cuIpcOpenEventHandle (CUevent* phEvent, CUipcEventHandle handle);

CUresult cuIpcGetMemHandle (CUipcMemHandle* pHandle, CUdeviceptr dptr);

CUresult cuIpcOpenMemHandle (CUdeviceptr* pdptr, CUipcMemHandle handle, uint Flags);

CUresult cuIpcCloseMemHandle (CUdeviceptr dptr);

CUresult cuMemHostRegister_v2 (void* p, size_t bytesize, uint Flags);

CUresult cuMemHostUnregister (void* p);

CUresult cuMemcpy (CUdeviceptr dst, CUdeviceptr src, size_t ByteCount);

CUresult cuMemcpyPeer (CUdeviceptr dstDevice, CUcontext dstContext, CUdeviceptr srcDevice, CUcontext srcContext, size_t ByteCount);

CUresult cuMemcpyHtoD_v2 (CUdeviceptr dstDevice, const(void)* srcHost, size_t ByteCount);

CUresult cuMemcpyDtoH_v2 (void* dstHost, CUdeviceptr srcDevice, size_t ByteCount);

CUresult cuMemcpyDtoD_v2 (CUdeviceptr dstDevice, CUdeviceptr srcDevice, size_t ByteCount);

CUresult cuMemcpyDtoA_v2 (CUarray dstArray, size_t dstOffset, CUdeviceptr srcDevice, size_t ByteCount);

CUresult cuMemcpyAtoD_v2 (CUdeviceptr dstDevice, CUarray srcArray, size_t srcOffset, size_t ByteCount);

CUresult cuMemcpyHtoA_v2 (CUarray dstArray, size_t dstOffset, const(void)* srcHost, size_t ByteCount);

CUresult cuMemcpyAtoH_v2 (void* dstHost, CUarray srcArray, size_t srcOffset, size_t ByteCount);

CUresult cuMemcpyAtoA_v2 (CUarray dstArray, size_t dstOffset, CUarray srcArray, size_t srcOffset, size_t ByteCount);

CUresult cuMemcpy2D_v2 (const(CUDA_MEMCPY2D)* pCopy);


CUresult cuMemcpy2DUnaligned_v2 (const(CUDA_MEMCPY2D)* pCopy);


CUresult cuMemcpy3D_v2 (const(CUDA_MEMCPY3D)* pCopy);

CUresult cuMemcpy3DPeer (const(CUDA_MEMCPY3D_PEER)* pCopy);

CUresult cuMemcpyAsync (CUdeviceptr dst, CUdeviceptr src, size_t ByteCount, CUstream hStream);

CUresult cuMemcpyPeerAsync (CUdeviceptr dstDevice, CUcontext dstContext, CUdeviceptr srcDevice, CUcontext srcContext, size_t ByteCount, CUstream hStream);

CUresult cuMemcpyHtoDAsync_v2 (CUdeviceptr dstDevice, const(void)* srcHost, size_t ByteCount, CUstream hStream);

CUresult cuMemcpyDtoHAsync_v2 (void* dstHost, CUdeviceptr srcDevice, size_t ByteCount, CUstream hStream);

CUresult cuMemcpyDtoDAsync_v2 (CUdeviceptr dstDevice, CUdeviceptr srcDevice, size_t ByteCount, CUstream hStream);

CUresult cuMemcpyHtoAAsync_v2 (CUarray dstArray, size_t dstOffset, const(void)* srcHost, size_t ByteCount, CUstream hStream);

CUresult cuMemcpyAtoHAsync_v2 (void* dstHost, CUarray srcArray, size_t srcOffset, size_t ByteCount, CUstream hStream);

CUresult cuMemcpy2DAsync_v2 (const(CUDA_MEMCPY2D)* pCopy, CUstream hStream);

CUresult cuMemcpy3DAsync_v2 (const(CUDA_MEMCPY3D)* pCopy, CUstream hStream);

CUresult cuMemcpy3DPeerAsync (const(CUDA_MEMCPY3D_PEER)* pCopy, CUstream hStream);

CUresult cuMemsetD8_v2 (CUdeviceptr dstDevice, ubyte uc, size_t N);

CUresult cuMemsetD16_v2 (CUdeviceptr dstDevice, ushort us, size_t N);

CUresult cuMemsetD32_v2 (CUdeviceptr dstDevice, uint ui, size_t N);

CUresult cuMemsetD2D8_v2 (CUdeviceptr dstDevice, size_t dstPitch, ubyte uc, size_t Width, size_t Height);

CUresult cuMemsetD2D16_v2 (CUdeviceptr dstDevice, size_t dstPitch, ushort us, size_t Width, size_t Height);

CUresult cuMemsetD2D32_v2 (CUdeviceptr dstDevice, size_t dstPitch, uint ui, size_t Width, size_t Height);

CUresult cuMemsetD8Async (CUdeviceptr dstDevice, ubyte uc, size_t N, CUstream hStream);

CUresult cuMemsetD16Async (CUdeviceptr dstDevice, ushort us, size_t N, CUstream hStream);

CUresult cuMemsetD32Async (CUdeviceptr dstDevice, uint ui, size_t N, CUstream hStream);

CUresult cuMemsetD2D8Async (CUdeviceptr dstDevice, size_t dstPitch, ubyte uc, size_t Width, size_t Height, CUstream hStream);

CUresult cuMemsetD2D16Async (CUdeviceptr dstDevice, size_t dstPitch, ushort us, size_t Width, size_t Height, CUstream hStream);

CUresult cuMemsetD2D32Async (CUdeviceptr dstDevice, size_t dstPitch, uint ui, size_t Width, size_t Height, CUstream hStream);

CUresult cuArrayCreate_v2 (CUarray* pHandle, const(CUDA_ARRAY_DESCRIPTOR)* pAllocateArray);

CUresult cuArrayGetDescriptor_v2 (CUDA_ARRAY_DESCRIPTOR* pArrayDescriptor, CUarray hArray);

CUresult cuArrayDestroy (CUarray hArray);

CUresult cuArray3DCreate_v2 (CUarray* pHandle, const(CUDA_ARRAY3D_DESCRIPTOR)* pAllocateArray);

CUresult cuArray3DGetDescriptor_v2 (CUDA_ARRAY3D_DESCRIPTOR* pArrayDescriptor, CUarray hArray);

CUresult cuMipmappedArrayCreate (CUmipmappedArray* pHandle, const(CUDA_ARRAY3D_DESCRIPTOR)* pMipmappedArrayDesc, uint numMipmapLevels);

CUresult cuMipmappedArrayGetLevel (CUarray* pLevelArray, CUmipmappedArray hMipmappedArray, uint level);

CUresult cuMipmappedArrayDestroy (CUmipmappedArray hMipmappedArray);

CUresult cuPointerGetAttribute (void* data, CUpointer_attribute attribute, CUdeviceptr ptr);

CUresult cuPointerSetAttribute (const(void)* value, CUpointer_attribute attribute, CUdeviceptr ptr);

CUresult cuPointerGetAttributes (uint numAttributes, CUpointer_attribute* attributes, void** data, CUdeviceptr ptr);

CUresult cuStreamCreate (CUstream* phStream, uint Flags);

CUresult cuStreamCreateWithPriority (CUstream* phStream, uint flags, int priority);

CUresult cuStreamGetPriority (CUstream hStream, int* priority);

CUresult cuStreamGetFlags (CUstream hStream, uint* flags);

CUresult cuStreamWaitEvent (CUstream hStream, CUevent hEvent, uint Flags);

CUresult cuStreamAddCallback (CUstream hStream, CUstreamCallback callback, void* userData, uint flags);

CUresult cuStreamAttachMemAsync (CUstream hStream, CUdeviceptr dptr, size_t length, uint flags);

CUresult cuStreamQuery (CUstream hStream);

CUresult cuStreamSynchronize (CUstream hStream);

CUresult cuStreamDestroy_v2 (CUstream hStream);

CUresult cuEventCreate (CUevent* phEvent, uint Flags);

CUresult cuEventRecord (CUevent hEvent, CUstream hStream);

CUresult cuEventQuery (CUevent hEvent);

CUresult cuEventSynchronize (CUevent hEvent);

CUresult cuEventDestroy_v2 (CUevent hEvent);

CUresult cuEventElapsedTime (float* pMilliseconds, CUevent hStart, CUevent hEnd);

CUresult cuFuncGetAttribute (int* pi, CUfunction_attribute attrib, CUfunction hfunc);

CUresult cuFuncSetCacheConfig (CUfunction hfunc, CUfunc_cache config);

CUresult cuFuncSetSharedMemConfig (CUfunction hfunc, CUsharedconfig config);

CUresult cuLaunchKernel (
    CUfunction f,
    uint gridDimX,
    uint gridDimY,
    uint gridDimZ,
    uint blockDimX,
    uint blockDimY,
    uint blockDimZ,
    uint sharedMemBytes,
    CUstream hStream,
    void** kernelParams,
    void** extra);

CUresult cuFuncSetBlockShape (CUfunction hfunc, int x, int y, int z);

CUresult cuFuncSetSharedSize (CUfunction hfunc, uint bytes);

CUresult cuParamSetSize (CUfunction hfunc, uint numbytes);

CUresult cuParamSeti (CUfunction hfunc, int offset, uint value);

CUresult cuParamSetf (CUfunction hfunc, int offset, float value);

CUresult cuParamSetv (CUfunction hfunc, int offset, void* ptr, uint numbytes);

CUresult cuLaunch (CUfunction f);

CUresult cuLaunchGrid (CUfunction f, int grid_width, int grid_height);

CUresult cuLaunchGridAsync (CUfunction f, int grid_width, int grid_height, CUstream hStream);

CUresult cuParamSetTexRef (CUfunction hfunc, int texunit, CUtexref hTexRef);

CUresult cuOccupancyMaxActiveBlocksPerMultiprocessor (int* numBlocks, CUfunction func, int blockSize, size_t dynamicSMemSize);

CUresult cuOccupancyMaxActiveBlocksPerMultiprocessorWithFlags (int* numBlocks, CUfunction func, int blockSize, size_t dynamicSMemSize, uint flags);

CUresult cuOccupancyMaxPotentialBlockSize (int* minGridSize, int* blockSize, CUfunction func, CUoccupancyB2DSize blockSizeToDynamicSMemSize, size_t dynamicSMemSize, int blockSizeLimit);

CUresult cuOccupancyMaxPotentialBlockSizeWithFlags (int* minGridSize, int* blockSize, CUfunction func, CUoccupancyB2DSize blockSizeToDynamicSMemSize, size_t dynamicSMemSize, int blockSizeLimit, uint flags);

CUresult cuTexRefSetArray (CUtexref hTexRef, CUarray hArray, uint Flags);

CUresult cuTexRefSetMipmappedArray (CUtexref hTexRef, CUmipmappedArray hMipmappedArray, uint Flags);

CUresult cuTexRefSetAddress_v2 (size_t* ByteOffset, CUtexref hTexRef, CUdeviceptr dptr, size_t bytes);

CUresult cuTexRefSetAddress2D_v3 (CUtexref hTexRef, const(CUDA_ARRAY_DESCRIPTOR)* desc, CUdeviceptr dptr, size_t Pitch);

CUresult cuTexRefSetFormat (CUtexref hTexRef, CUarray_format fmt, int NumPackedComponents);

CUresult cuTexRefSetAddressMode (CUtexref hTexRef, int dim, CUaddress_mode am);

CUresult cuTexRefSetFilterMode (CUtexref hTexRef, CUfilter_mode fm);

CUresult cuTexRefSetMipmapFilterMode (CUtexref hTexRef, CUfilter_mode fm);

CUresult cuTexRefSetMipmapLevelBias (CUtexref hTexRef, float bias);

CUresult cuTexRefSetMipmapLevelClamp (CUtexref hTexRef, float minMipmapLevelClamp, float maxMipmapLevelClamp);

CUresult cuTexRefSetMaxAnisotropy (CUtexref hTexRef, uint maxAniso);

CUresult cuTexRefSetFlags (CUtexref hTexRef, uint Flags);

CUresult cuTexRefGetAddress_v2 (CUdeviceptr* pdptr, CUtexref hTexRef);

CUresult cuTexRefGetArray (CUarray* phArray, CUtexref hTexRef);

CUresult cuTexRefGetMipmappedArray (CUmipmappedArray* phMipmappedArray, CUtexref hTexRef);

CUresult cuTexRefGetAddressMode (CUaddress_mode* pam, CUtexref hTexRef, int dim);

CUresult cuTexRefGetFilterMode (CUfilter_mode* pfm, CUtexref hTexRef);

CUresult cuTexRefGetFormat (CUarray_format* pFormat, int* pNumChannels, CUtexref hTexRef);

CUresult cuTexRefGetMipmapFilterMode (CUfilter_mode* pfm, CUtexref hTexRef);

CUresult cuTexRefGetMipmapLevelBias (float* pbias, CUtexref hTexRef);

CUresult cuTexRefGetMipmapLevelClamp (float* pminMipmapLevelClamp, float* pmaxMipmapLevelClamp, CUtexref hTexRef);

CUresult cuTexRefGetMaxAnisotropy (int* pmaxAniso, CUtexref hTexRef);

CUresult cuTexRefGetFlags (uint* pFlags, CUtexref hTexRef);

CUresult cuTexRefCreate (CUtexref* pTexRef);

CUresult cuTexRefDestroy (CUtexref hTexRef);

CUresult cuSurfRefSetArray (CUsurfref hSurfRef, CUarray hArray, uint Flags);

CUresult cuSurfRefGetArray (CUarray* phArray, CUsurfref hSurfRef);

CUresult cuTexObjectCreate (CUtexObject* pTexObject, const(CUDA_RESOURCE_DESC)* pResDesc, const(CUDA_TEXTURE_DESC)* pTexDesc, const(CUDA_RESOURCE_VIEW_DESC)* pResViewDesc);

CUresult cuTexObjectDestroy (CUtexObject texObject);

CUresult cuTexObjectGetResourceDesc (CUDA_RESOURCE_DESC* pResDesc, CUtexObject texObject);

CUresult cuTexObjectGetTextureDesc (CUDA_TEXTURE_DESC* pTexDesc, CUtexObject texObject);

CUresult cuTexObjectGetResourceViewDesc (CUDA_RESOURCE_VIEW_DESC* pResViewDesc, CUtexObject texObject);

CUresult cuSurfObjectCreate (CUsurfObject* pSurfObject, const(CUDA_RESOURCE_DESC)* pResDesc);

CUresult cuSurfObjectDestroy (CUsurfObject surfObject);

CUresult cuSurfObjectGetResourceDesc (CUDA_RESOURCE_DESC* pResDesc, CUsurfObject surfObject);

CUresult cuDeviceCanAccessPeer (int* canAccessPeer, CUdevice dev, CUdevice peerDev);

CUresult cuCtxEnablePeerAccess (CUcontext peerContext, uint Flags);

CUresult cuCtxDisablePeerAccess (CUcontext peerContext);

CUresult cuGraphicsUnregisterResource (CUgraphicsResource resource);

CUresult cuGraphicsSubResourceGetMappedArray (CUarray* pArray, CUgraphicsResource resource, uint arrayIndex, uint mipLevel);

CUresult cuGraphicsResourceGetMappedMipmappedArray (CUmipmappedArray* pMipmappedArray, CUgraphicsResource resource);

CUresult cuGraphicsResourceGetMappedPointer_v2 (CUdeviceptr* pDevPtr, size_t* pSize, CUgraphicsResource resource);

CUresult cuGraphicsResourceSetMapFlags_v2 (CUgraphicsResource resource, uint flags);


CUresult cuGraphicsMapResources (uint count, CUgraphicsResource* resources, CUstream hStream);

CUresult cuGraphicsUnmapResources (uint count, CUgraphicsResource* resources, CUstream hStream);

/** @} */ /* END CUDA_GRAPHICS */

CUresult cuGetExportTable (const(void*)* ppExportTable, const(CUuuid)* pExportTableId);


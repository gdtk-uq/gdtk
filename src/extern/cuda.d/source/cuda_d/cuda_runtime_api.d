/*
  D bindings for CUDA.
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module cuda_d.cuda_runtime_api;

import cuda_d.cuda;
import cuda_d.cublas_api;
import cuda_d.vector_types;

extern (C):

enum CUDART_VERSION = 7050;

extern (D) auto __CUDART_API_PTDS(T)(auto ref T api)
{
    return api;
}

extern (D) auto __CUDART_API_PTSZ(T)(auto ref T api)
{
    return api;
}

cudaError_t cudaDeviceReset ();

cudaError_t cudaDeviceSynchronize ();

cudaError_t cudaDeviceSetLimit (cudaLimit limit, size_t value);

cudaError_t cudaDeviceGetLimit (size_t* pValue, cudaLimit limit);

cudaError_t cudaDeviceGetCacheConfig (cudaFuncCache* pCacheConfig);

cudaError_t cudaDeviceGetStreamPriorityRange (int* leastPriority, int* greatestPriority);

cudaError_t cudaDeviceSetCacheConfig (cudaFuncCache cacheConfig);

cudaError_t cudaDeviceGetSharedMemConfig (cudaSharedMemConfig* pConfig);

cudaError_t cudaDeviceSetSharedMemConfig (cudaSharedMemConfig config);

cudaError_t cudaDeviceGetByPCIBusId (int* device, const(char)* pciBusId);

cudaError_t cudaDeviceGetPCIBusId (char* pciBusId, int len, int device);

cudaError_t cudaIpcGetEventHandle (cudaIpcEventHandle_t* handle, cudaEvent_t event);

cudaError_t cudaIpcOpenEventHandle (cudaEvent_t* event, cudaIpcEventHandle_t handle);

cudaError_t cudaIpcGetMemHandle (cudaIpcMemHandle_t* handle, void* devPtr);

cudaError_t cudaIpcOpenMemHandle (void** devPtr, cudaIpcMemHandle_t handle, uint flags);

cudaError_t cudaIpcCloseMemHandle (void* devPtr);

cudaError_t cudaThreadExit ();

cudaError_t cudaThreadSynchronize ();

cudaError_t cudaThreadSetLimit (cudaLimit limit, size_t value);

cudaError_t cudaThreadGetLimit (size_t* pValue, cudaLimit limit);

cudaError_t cudaThreadGetCacheConfig (cudaFuncCache* pCacheConfig);

cudaError_t cudaThreadSetCacheConfig (cudaFuncCache cacheConfig);

cudaError_t cudaGetLastError ();

cudaError_t cudaPeekAtLastError ();

const(char)* cudaGetErrorName (cudaError_t error);


const(char)* cudaGetErrorString (cudaError_t error);

cudaError_t cudaGetDeviceCount (int* count);

cudaError_t cudaGetDeviceProperties (cudaDeviceProp* prop, int device);

//cudaError_t cudaDeviceGetAttribute (int* value, cudaDeviceAttr attr, int device);

cudaError_t cudaChooseDevice (int* device, const(cudaDeviceProp)* prop);

cudaError_t cudaSetDevice (int device);

cudaError_t cudaGetDevice (int* device);

cudaError_t cudaSetValidDevices (int* device_arr, int len);

cudaError_t cudaSetDeviceFlags (uint flags);

cudaError_t cudaGetDeviceFlags (uint* flags);

cudaError_t cudaStreamCreate (cudaStream_t* pStream);

cudaError_t cudaStreamCreateWithFlags (cudaStream_t* pStream, uint flags);

cudaError_t cudaStreamCreateWithPriority (cudaStream_t* pStream, uint flags, int priority);

cudaError_t cudaStreamGetPriority (cudaStream_t hStream, int* priority);

cudaError_t cudaStreamGetFlags (cudaStream_t hStream, uint* flags);

cudaError_t cudaStreamDestroy (cudaStream_t stream);

cudaError_t cudaStreamWaitEvent (cudaStream_t stream, cudaEvent_t event, uint flags);

alias cudaStreamCallback_t = void function (cudaStream_t stream, cudaError_t status, void* userData);


cudaError_t cudaStreamAddCallback (
    cudaStream_t stream,
    cudaStreamCallback_t callback,
    void* userData,
    uint flags);


cudaError_t cudaStreamSynchronize (cudaStream_t stream);

cudaError_t cudaStreamQuery (cudaStream_t stream);

cudaError_t cudaStreamAttachMemAsync (cudaStream_t stream, void* devPtr, size_t length, uint flags);

cudaError_t cudaEventCreate (cudaEvent_t* event);

cudaError_t cudaEventCreateWithFlags (cudaEvent_t* event, uint flags);

cudaError_t cudaEventRecord (cudaEvent_t event, cudaStream_t stream);


cudaError_t cudaEventQuery (cudaEvent_t event);

cudaError_t cudaEventSynchronize (cudaEvent_t event);

cudaError_t cudaEventDestroy (cudaEvent_t event);

cudaError_t cudaEventElapsedTime (float* ms, cudaEvent_t start, cudaEvent_t end);

cudaError_t cudaLaunchKernel (const(void)* func, dim3 gridDim, dim3 blockDim, void** args, size_t sharedMem, cudaStream_t stream);

cudaError_t cudaFuncSetCacheConfig (const(void)* func, cudaFuncCache cacheConfig);

cudaError_t cudaFuncSetSharedMemConfig (const(void)* func, cudaSharedMemConfig config);

//cudaError_t cudaFuncGetAttributes (cudaFuncAttributes* attr, const(void)* func);

cudaError_t cudaSetDoubleForDevice (double* d);

cudaError_t cudaSetDoubleForHost (double* d);

cudaError_t cudaOccupancyMaxActiveBlocksPerMultiprocessor (int* numBlocks, const(void)* func, int blockSize, size_t dynamicSMemSize);

cudaError_t cudaOccupancyMaxActiveBlocksPerMultiprocessorWithFlags (int* numBlocks, const(void)* func, int blockSize, size_t dynamicSMemSize, uint flags);

cudaError_t cudaConfigureCall (dim3 gridDim, dim3 blockDim, size_t sharedMem, cudaStream_t stream);

cudaError_t cudaSetupArgument (const(void)* arg, size_t size, size_t offset);

cudaError_t cudaLaunch (const(void)* func);

cudaError_t cudaMallocManaged (void** devPtr, size_t size, uint flags);

cudaError_t cudaMalloc (void** devPtr, size_t size);

cudaError_t cudaMallocHost (void** ptr, size_t size);

cudaError_t cudaMallocPitch (void** devPtr, size_t* pitch, size_t width, size_t height);

//cudaError_t cudaMallocArray (cudaArray_t* array, const(cudaChannelFormatDesc)* desc, size_t width, size_t height, uint flags);

cudaError_t cudaFree (void* devPtr);

cudaError_t cudaFreeHost (void* ptr);

cudaError_t cudaFreeArray (cudaArray_t array);

//cudaError_t cudaFreeMipmappedArray (cudaMipmappedArray_t mipmappedArray);

cudaError_t cudaHostAlloc (void** pHost, size_t size, uint flags);

cudaError_t cudaHostRegister (void* ptr, size_t size, uint flags);

cudaError_t cudaHostUnregister (void* ptr);

cudaError_t cudaHostGetDevicePointer (void** pDevice, void* pHost, uint flags);

cudaError_t cudaHostGetFlags (uint* pFlags, void* pHost);

//cudaError_t cudaMalloc3D (cudaPitchedPtr* pitchedDevPtr, cudaExtent extent);

//cudaError_t cudaMalloc3DArray (cudaArray_t* array, const(cudaChannelFormatDesc)* desc, cudaExtent extent, uint flags);

//cudaError_t cudaMallocMipmappedArray (cudaMipmappedArray_t* mipmappedArray, const(cudaChannelFormatDesc)* desc, cudaExtent extent, uint numLevels, uint flags);

//cudaError_t cudaGetMipmappedArrayLevel (cudaArray_t* levelArray, cudaMipmappedArray_const_t mipmappedArray, uint level);

//http://horacio9573.no-ip.org/cuda/structcudaMemcpy3DParms.html
//cudaError_t cudaMemcpy3D (const(cudaMemcpy3DParms)* p);

//cudaError_t cudaMemcpy3DPeer (const(cudaMemcpy3DPeerParms)* p);

//cudaError_t cudaMemcpy3DAsync (const(cudaMemcpy3DParms)* p, cudaStream_t stream);

//cudaError_t cudaMemcpy3DPeerAsync (const(cudaMemcpy3DPeerParms)* p, cudaStream_t stream);

cudaError_t cudaMemGetInfo (size_t* free, size_t* total);

//cudaError_t cudaArrayGetInfo (cudaChannelFormatDesc* desc, cudaExtent* extent, uint* flags, cudaArray_t array);

cudaError_t cudaMemcpy (void* dst, const(void)* src, size_t count, cudaMemcpyKind kind);

cudaError_t cudaMemcpyPeer (void* dst, int dstDevice, const(void)* src, int srcDevice, size_t count);

cudaError_t cudaMemcpyToArray (cudaArray_t dst, size_t wOffset, size_t hOffset, const(void)* src, size_t count, cudaMemcpyKind kind);

cudaError_t cudaMemcpyFromArray (void* dst, cudaArray_const_t src, size_t wOffset, size_t hOffset, size_t count, cudaMemcpyKind kind);

cudaError_t cudaMemcpyArrayToArray (cudaArray_t dst, size_t wOffsetDst, size_t hOffsetDst, cudaArray_const_t src, size_t wOffsetSrc, size_t hOffsetSrc, size_t count, cudaMemcpyKind kind);

cudaError_t cudaMemcpy2D (void* dst, size_t dpitch, const(void)* src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind);

cudaError_t cudaMemcpy2DToArray (cudaArray_t dst, size_t wOffset, size_t hOffset, const(void)* src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind);

cudaError_t cudaMemcpy2DFromArray (void* dst, size_t dpitch, cudaArray_const_t src, size_t wOffset, size_t hOffset, size_t width, size_t height, cudaMemcpyKind kind);

cudaError_t cudaMemcpy2DArrayToArray (cudaArray_t dst, size_t wOffsetDst, size_t hOffsetDst, cudaArray_const_t src, size_t wOffsetSrc, size_t hOffsetSrc, size_t width, size_t height, cudaMemcpyKind kind);

cudaError_t cudaMemcpyToSymbol (const(void)* symbol, const(void)* src, size_t count, size_t offset, cudaMemcpyKind kind);

cudaError_t cudaMemcpyFromSymbol (void* dst, const(void)* symbol, size_t count, size_t offset, cudaMemcpyKind kind);

cudaError_t cudaMemcpyAsync (void* dst, const(void)* src, size_t count, cudaMemcpyKind kind, cudaStream_t stream);

cudaError_t cudaMemcpyPeerAsync (void* dst, int dstDevice, const(void)* src, int srcDevice, size_t count, cudaStream_t stream);

cudaError_t cudaMemcpyToArrayAsync (cudaArray_t dst, size_t wOffset, size_t hOffset, const(void)* src, size_t count, cudaMemcpyKind kind, cudaStream_t stream);

cudaError_t cudaMemcpyFromArrayAsync (void* dst, cudaArray_const_t src, size_t wOffset, size_t hOffset, size_t count, cudaMemcpyKind kind, cudaStream_t stream);

cudaError_t cudaMemcpy2DAsync (void* dst, size_t dpitch, const(void)* src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind, cudaStream_t stream);

cudaError_t cudaMemcpy2DToArrayAsync (cudaArray_t dst, size_t wOffset, size_t hOffset, const(void)* src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind, cudaStream_t stream);

cudaError_t cudaMemcpy2DFromArrayAsync (void* dst, size_t dpitch, cudaArray_const_t src, size_t wOffset, size_t hOffset, size_t width, size_t height, cudaMemcpyKind kind, cudaStream_t stream);

cudaError_t cudaMemcpyToSymbolAsync (const(void)* symbol, const(void)* src, size_t count, size_t offset, cudaMemcpyKind kind, cudaStream_t stream);

cudaError_t cudaMemcpyFromSymbolAsync (void* dst, const(void)* symbol, size_t count, size_t offset, cudaMemcpyKind kind, cudaStream_t stream);

cudaError_t cudaMemset (void* devPtr, int value, size_t count);

cudaError_t cudaMemset2D (void* devPtr, size_t pitch, int value, size_t width, size_t height);

//cudaError_t cudaMemset3D (cudaPitchedPtr pitchedDevPtr, int value, cudaExtent extent);

//cudaError_t cudaMemsetAsync (void* devPtr, int value, size_t count, cudaStream_t stream);

//cudaError_t cudaMemset2DAsync (void* devPtr, size_t pitch, int value, size_t width, size_t height, cudaStream_t stream);

//cudaError_t cudaMemset3DAsync (cudaPitchedPtr pitchedDevPtr, int value, cudaExtent extent, cudaStream_t stream);

//cudaError_t cudaGetSymbolAddress (void** devPtr, const(void)* symbol);

//cudaError_t cudaGetSymbolSize (size_t* size, const(void)* symbol);

//cudaError_t cudaPointerGetAttributes (cudaPointerAttributes* attributes, const(void)* ptr);

//cudaError_t cudaDeviceCanAccessPeer (int* canAccessPeer, int device, int peerDevice);

//cudaError_t cudaDeviceEnablePeerAccess (int peerDevice, uint flags);

//cudaError_t cudaDeviceDisablePeerAccess (int peerDevice);

//cudaChannelFormatDesc cudaCreateChannelDesc (int x, int y, int z, int w, cudaChannelFormatKind f);

//cudaError_t cudaGraphicsUnregisterResource (cudaGraphicsResource_t resource);

//cudaError_t cudaGraphicsResourceSetMapFlags (cudaGraphicsResource_t resource, uint flags);

//cudaError_t cudaGraphicsMapResources (int count, cudaGraphicsResource_t* resources, cudaStream_t stream);

//cudaError_t cudaGraphicsUnmapResources (int count, cudaGraphicsResource_t* resources, cudaStream_t stream);

//cudaError_t cudaGraphicsResourceGetMappedPointer (void** devPtr, size_t* size, cudaGraphicsResource_t resource);

//cudaError_t cudaGraphicsSubResourceGetMappedArray (cudaArray_t* array, cudaGraphicsResource_t resource, uint arrayIndex, uint mipLevel);

//cudaError_t cudaGraphicsResourceGetMappedMipmappedArray (cudaMipmappedArray_t* mipmappedArray, cudaGraphicsResource_t resource);

//cudaError_t cudaGetChannelDesc (cudaChannelFormatDesc* desc, cudaArray_const_t array);

//cudaError_t cudaBindTexture (size_t* offset, const(textureReference)* texref, const(void)* devPtr, const(cudaChannelFormatDesc)* desc, size_t size);

//cudaError_t cudaBindTexture2D (size_t* offset, const(textureReference)* texref, const(void)* devPtr, const(cudaChannelFormatDesc)* desc, size_t width, size_t height, size_t pitch);

//cudaError_t cudaBindTextureToArray (const(textureReference)* texref, cudaArray_const_t array, const(cudaChannelFormatDesc)* desc);

//cudaError_t cudaBindTextureToMipmappedArray (const(textureReference)* texref, cudaMipmappedArray_const_t mipmappedArray, const(cudaChannelFormatDesc)* desc);

//cudaError_t cudaUnbindTexture (const(textureReference)* texref);

//cudaError_t cudaGetTextureAlignmentOffset (size_t* offset, const(textureReference)* texref);

//cudaError_t cudaGetTextureReference (const(textureReference*)* texref, const(void)* symbol);

//cudaError_t cudaBindSurfaceToArray (const(surfaceReference)* surfref, cudaArray_const_t array, const(cudaChannelFormatDesc)* desc);

//cudaError_t cudaGetSurfaceReference (const(surfaceReference*)* surfref, const(void)* symbol);

//cudaError_t cudaCreateTextureObject (cudaTextureObject_t* pTexObject, const(cudaResourceDesc)* pResDesc, const(cudaTextureDesc)* pTexDesc, const(cudaResourceViewDesc)* pResViewDesc);

//cudaError_t cudaDestroyTextureObject (cudaTextureObject_t texObject);

//cudaError_t cudaGetTextureObjectResourceDesc (cudaResourceDesc* pResDesc, cudaTextureObject_t texObject);

//cudaError_t cudaGetTextureObjectTextureDesc (cudaTextureDesc* pTexDesc, cudaTextureObject_t texObject);

//cudaError_t cudaGetTextureObjectResourceViewDesc (cudaResourceViewDesc* pResViewDesc, cudaTextureObject_t texObject);

//cudaError_t cudaCreateSurfaceObject (cudaSurfaceObject_t* pSurfObject, const(cudaResourceDesc)* pResDesc);

//cudaError_t cudaDestroySurfaceObject (cudaSurfaceObject_t surfObject);

//cudaError_t cudaGetSurfaceObjectResourceDesc (cudaResourceDesc* pResDesc, cudaSurfaceObject_t surfObject);

//cudaError_t cudaDriverGetVersion (int* driverVersion);

//cudaError_t cudaRuntimeGetVersion (int* runtimeVersion);

//cudaError_t cudaGetExportTable (const(void*)* ppExportTable, const(cudaUUID_t)* pExportTableId);

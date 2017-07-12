# makefile for the CUDA D Hello World test
# to run, execute: $LD_LIBRARY_PATH=$LD_LIBRARY_PATH:. ./hello

DMD ?= dmd
NVCC ?= nvcc
GPP ?= g++

CUDA_DIR := ../source/cuda_d
CUDA_FILES := $(CUDA_DIR)/cublas_api.d \
	$(CUDA_DIR)/cublasXt.d \
	$(CUDA_DIR)/cudalibxt.d \
	$(CUDA_DIR)/curand.d \
	$(CUDA_DIR)/cublas.d \
	$(CUDA_DIR)/cuComplex.d \
	$(CUDA_DIR)/cuda_profiler_api.d \
	$(CUDA_DIR)/vector_types.d \
	$(CUDA_DIR)/cublas_v2.d \
	$(CUDA_DIR)/cuda.d \
	$(CUDA_DIR)/cuda_runtime_api.d

PROGRAMS := hello

clean:
	- rm *.o
	- rm *.so
	- rm *~
	- rm $(PROGRAMS)

install: hello_world.d hello_world_kernel.cu
	$(NVCC) --shared -o libhwkernel.so hello_world_kernel.cu --compiler-options '-fPIC'  
	$(DMD) hello_world.d -I$(CUDA_DIR)/ $(CUDA_FILES) -L-L. -L-L/usr/local/cuda-7.5/lib64 \
		-L-lcuda -L-lcudart -L-lphobos2 -L-lhwkernel

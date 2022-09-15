# makefile for chicken
# PJ 2022-09-11

SOURCE_FILES := main.cu block.cu cell.cu face.cu vertex.cu flux.cu flow.cu gas.cu vector3.cu number.cu

# Default build is for double numbers.
chicken : $(SOURCE_FILES)
	nvcc -DCUDA -std c++17 -lineinfo -o chicken main.cu

chicken_f : $(SOURCE_FILES)
	nvcc -DCUDA -DFLOAT_NUMBERS -std c++17 -lineinfo -o chicken main.cu

chicken_host_only : $(SOURCE_FILES)
	nvcc -lineinfo -std c++17 -o chicken main.cu

chicken_run : chicken
	./chicken

# Note that profiling needs to be done via sudo
# because we need access to the hardware counters.
chicken_profile : chicken
	/usr/local/cuda-11.7/bin/ncu --set full --clock-control none --target-processes all \
		-o chicken -f ./chicken

clean :
	- rm chicken chicken.ncu-rep

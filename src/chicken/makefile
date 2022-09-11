# makefile for chicken
# PJ 2022-09-11

SOURCE_FILES := main.cu cell.cu face.cu flux.cu flow.cu gas.cu vector3.cu number.cu

# Default build is for double numbers.
chicken : $(SOURCE_FILES)
	nvcc -lineinfo -o chicken main.cu

chicken_f : $(SOURCE_FILES)
	nvcc -DFLOAT_NUMBERS -lineinfo -o chicken main.cu

chicken_run : chicken
	./chicken

# Note that profiling needs to be done via sudo
# because we need access to the hardware counters.
chicken_profile : chicken
	/usr/local/cuda-11.7/bin/ncu --set full --clock-control none --target-processes all \
		-o chicken -f ./chicken

clean :
	- rm chicken chicken.ncu-rep

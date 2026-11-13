# makefile for the tests
# PJ 2022-11-19 Finally got tired of typing these commands manually.
#
PROGRAMS := gas_test rsla_test vector3_test spline_test

test: $(PROGRAMS)
	./gas_test
	@echo "--------------------------------------"
	./rsla_test
	@echo "--------------------------------------"
	./vector3_test
	@echo "--------------------------------------"
	./spline_test
	@echo "-Done---------------------------------"

clean:
	-rm $(PROGRAMS)

ifeq ($(HIP),1)
# HIP (AMD GPU) build targets
gas_test: gas_test.cu ../number.cu ../gas.cu
	hipcc -std=c++17 -o gas_test -DIDEAL_AIR gas_test.cu

rsla_test: rsla_test.cu ../number.cu ../rsla.cu
	hipcc -std=c++17 -o rsla_test rsla_test.cu

vector3_test: vector3_test.cu ../vector3.cu ../number.cu
	hipcc -std=c++17 -o vector3_test vector3_test.cu

spline_test: spline_test.cu ../spline.cu ../rsla.cu ../number.cu
	hipcc -std=c++17 -o spline_test spline_test.cu
else
# CUDA (NVIDIA GPU) build targets
gas_test: gas_test.cu ../number.cu ../gas.cu
	nvcc --expt-relaxed-constexpr -o gas_test -DIDEAL_AIR gas_test.cu

rsla_test: rsla_test.cu ../number.cu ../rsla.cu
	nvcc --expt-relaxed-constexpr -o rsla_test rsla_test.cu

vector3_test: vector3_test.cu ../vector3.cu ../number.cu
	nvcc --expt-relaxed-constexpr -o vector3_test vector3_test.cu

spline_test: spline_test.cu ../spline.cu ../rsla.cu ../number.cu
	nvcc --expt-relaxed-constexpr -o spline_test spline_test.cu
endif

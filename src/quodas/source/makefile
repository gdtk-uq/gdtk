# makefile: quodas
# Kyle A. Damm (27-01-2020)

DMD := ldmd2
VERSION ?= fast
ifeq ($(VERSION), debug)
	DFLAGS := -w -g -debug
        OF := -of
endif
ifeq ($(VERSION), fast)
	DFLAGS := -w -g -O -release -inline -boundscheck=off
        OF := -of
endif
WITH_COMPLEX_NUMBERS ?= 0
ifeq ($(WITH_COMPLEX_NUMBERS), 1)
	CMPLXFLAG := -version=complex_numbers
else
	CMPLXFLAG := -version=real_numbers
endif

FILES := main.d \
	complex.d \
	number.d \
	config.d \
	block.d \
	finite_volume.d \
	flux.d \
	flow_solver.d \
	optimizer.d \
	linalg.d \
	derivative.d

quodas: $(FILES)
	$(DMD) $(DFLAGS) $(OF)$@ $(FILES) $(CMPLXFLAG)

complex_number_test: complexify.d number.d
	$(DMD) $(DFLAGS) $(OF)$@ complexify.d number.d -version=complex_numbers -version=complex_number_test

matrix_test: matrix.d number.d complexify.d
	$(DMD) $(DFLAGS) $(OF)$@ matrix.d number.d complexify.d -version=matrix_test

clean:
	- rm *.o
	- rm test
	- rm quodas
	- rm complex_number_test
	- rm matrix_test

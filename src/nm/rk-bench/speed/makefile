# makefile for the Runge-Kutta benchmark programs (from 2022)
# PJ 2025-09-01 extracted from the nm/makefile
#

bench: rungekutta_bench_d rungekutta_bench.go rungekutta_bench_cr
	@echo "Start benchmark programs"
	./rungekutta_bench_d
	go run rungekutta_bench.go
	./rungekutta_bench_cr
	@echo "Finished benchmark programs"

clean:
	- rm -f *.o *.obj
	- rm -f ./rungekutta_bench_d ./rungekutta_bench_cr

# Specific targets, for individual artefacts.

rungekutta_bench_d: rungekutta_bench.d
	ldc2 -of$@ -g -O -boundscheck=off rungekutta_bench.d

rungekutta_bench_cr: rungekutta_bench.cr
	crystal build -o$@ --release rungekutta_bench.cr


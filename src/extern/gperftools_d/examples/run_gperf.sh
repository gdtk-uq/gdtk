dub --compiler=ldc2
CPUPROFILE=./prof.out ./profiler_example
pprof --gv ./profiler_example ./prof.out
# pprof --pdf ./profiler_example ./prof.out > profile.pdf
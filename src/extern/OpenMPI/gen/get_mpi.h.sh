set -e

echo '#include <mpi.h>' > test.c
echo $(mpicc -M test.c | grep -o -E '[^[:space:]]+mpi.h')
rm test.*

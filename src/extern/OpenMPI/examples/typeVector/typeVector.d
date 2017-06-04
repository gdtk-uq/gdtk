///Found on http://geco.mines.edu/workshop/class2/examples/mpi/index.html
//TODO typeVector check to ensure it works properly...
/*
Shows how to use MPI_Type_vector to send noncontiguous blocks of data
 and MPI_Get_count and MPI_Get_elements to see the number of elements sent
*/
import core.stdc.stdio;
import core.stdc.stdlib;
import std.algorithm;
import std.array;
import std.string;
import core.memory;
import mpi;
import mpi.util;

int main(string[] args)
{
    int argc = cast(int)args.length;
    auto argv = args.toArgv();

    int myid, numprocs, mpi_err;
    immutable SIZE = 25;
    double[SIZE] svect;
    double[SIZE] rvect;
    int i, bonk1, bonk2, numx, stride, extent;
    MPI_Datatype MPI_LEFT_RITE;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    stride=5;
    numx=(SIZE+1)/stride;
    extent=1;
    if(myid == 1)
    {
      printf("numx=%d  extent=%d stride=%d\n", numx, extent, stride);
    }
    mpi_err=MPI_Type_vector(numx, extent, stride, MPI_DOUBLE, &MPI_LEFT_RITE);
    mpi_err=MPI_Type_commit(&MPI_LEFT_RITE);
    if(myid == 0)
    {
        for (i=0; i<SIZE; i++)
            svect[i]=i;
        MPI_Send(svect.ptr, 1, MPI_LEFT_RITE, 1, 100, MPI_COMM_WORLD);
    }
    if(myid == 1)
    {
        for (i=0; i<SIZE; i++)
            rvect[i]=-1;
        MPI_Recv(rvect.ptr, 1, MPI_LEFT_RITE, 0, 100, MPI_COMM_WORLD, &status);
    }
    if(myid == 1)
    {
        MPI_Get_count(&status, MPI_LEFT_RITE, &bonk1);
        MPI_Get_elements(&status, MPI_DOUBLE, &bonk2);
        printf("got %d elements of type MY_TYPE\n", bonk1);
        printf("which contained %d elements of type MPI_DOUBLE\n", bonk2);
        for (i=0; i<SIZE; i++)
            if(rvect[i] != -1)
                printf("%d %g\n", i, rvect[i]);
    }
    return MPI_Finalize();
}
/*
output
numx=5  extent=5 stride=1
got 1 elements of type MY_TYPE
which contained 5 elements of type MPI_DOUBLE
0 0
5 5
10 10
15 15
20 20
*/

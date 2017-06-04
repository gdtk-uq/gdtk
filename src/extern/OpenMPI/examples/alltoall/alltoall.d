///Found on http://geco.mines.edu/workshop/class2/examples/mpi/index.html

/*
! This program shows how to use MPI_Alltoall.  Each processor
! send/rec a different  random number to/from other processors.
*/
import core.stdc.stdio;
import std.algorithm;
import std.array;
import std.string;
import std.random;
import core.memory;
import mpi;
import mpi.util;

/* globals */
int numnodes, myid, mpi_err;
immutable mpi_root = 0;
/* end module  */

void init_it(int* argc, char*** argv)
{
    mpi_err = MPI_Init(argc, argv);
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &numnodes );
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}

int main(string[] args)
{
    int argc = cast(int)args.length;
    auto argv = args.toArgv();

    int* sray, rray;
    int* sdisp, scounts, rdisp, rcounts;
    int ssize, rsize, i, k, j;
    float z;

    init_it(&argc, &argv);
    scounts= cast(int*)GC.malloc(int.sizeof*numnodes, GC.BlkAttr.NO_SCAN);
    rcounts= cast(int*)GC.malloc(int.sizeof*numnodes, GC.BlkAttr.NO_SCAN);
    sdisp=   cast(int*)GC.malloc(int.sizeof*numnodes, GC.BlkAttr.NO_SCAN);
    rdisp=   cast(int*)GC.malloc(int.sizeof*numnodes, GC.BlkAttr.NO_SCAN);
/*
! seed the random number generator with a
! different number on each processor
*/
//    seed_random(myid); not needed
/* find  data to send */
    for(i=0; i<numnodes; i++)
    {
        scounts[i] = uniform(0,100);
    }
    printf("myid= %d scounts=", myid);
    for(i=0; i<numnodes; i++)
        printf("%d ",scounts[i]);
    printf("\n");
/* send the data */
    mpi_err = MPI_Alltoall( scounts,1,MPI_INT,
                            rcounts,1,MPI_INT,
                            MPI_COMM_WORLD);
    printf("myid= %d rcounts=",myid);
    for(i=0; i<numnodes; i++)
        printf("%d ",rcounts[i]);
    printf("\n");
    return MPI_Finalize();
}


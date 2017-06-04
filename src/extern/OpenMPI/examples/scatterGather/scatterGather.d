///Found on http://geco.mines.edu/workshop/class2/examples/mpi/index.html
/*
! This program shows how to use MPI_Scatter and MPI_Gather
! Each processor gets different data from the root processor
! by way of mpi_scatter.  The data is summed and then sent back
! to the root processor using MPI_Gather.  The root processor
! then prints the global sum.
*/
import core.stdc.stdio;
import std.algorithm;
import std.array;
import std.string;
import core.memory;
import mpi;
import mpi.util;

/*  globals */
int numnodes, myid, mpi_err;
immutable mpi_root = 0;
/* end globals  */

void init_it(int *argc, char ***argv)
{
    mpi_err = MPI_Init(argc,argv);
    mpi_err = MPI_Comm_size( MPI_COMM_WORLD, &numnodes );
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}

int main(string[] args)
{
    int argc = cast(int)args.length;
    auto argv = args.toArgv();

    int* myray, send_ray, back_ray;
    int count;
    int size, mysize, i, k, j, total;

    init_it(&argc, &argv);
/* each processor will get count elements from the root */
    count = 4;
    myray = cast(int*)GC.malloc(count*int.sizeof, GC.BlkAttr.NO_SCAN);
/* create the data to be sent on the root */
    if(myid == mpi_root)
    {
        size=count*numnodes;
        send_ray = cast(int*)GC.malloc(size*int.sizeof, GC.BlkAttr.NO_SCAN);
        back_ray = cast(int*)GC.malloc(numnodes*int.sizeof, GC.BlkAttr.NO_SCAN);
        for(i=0; i<size; i++)
            send_ray[i] = i;
    }
/* send different data to each processor */
    mpi_err = MPI_Scatter(  send_ray, count,   MPI_INT,
                            myray,    count,   MPI_INT,
                            mpi_root,
                            MPI_COMM_WORLD);

/* each processor does a local sum */
    total=0;
    for(i=0; i<count; i++)
        total = total+myray[i];
    printf("myid= %d total= %d\n", myid, total);
/* send the local sums back to the root */
    mpi_err = MPI_Gather(&total,  1,  MPI_INT,
                        back_ray, 1,  MPI_INT,
                        mpi_root,
                        MPI_COMM_WORLD);
/* the root prints the global sum */
    if(myid == mpi_root)
    {
      total=0;
      for(i=0; i<numnodes; i++)
        total = total + back_ray[i];
      printf("results from all processors= %d\n", total);
    }
    return MPI_Finalize();
}

///Found on http://geco.mines.edu/workshop/class2/examples/mpi/index.html

import core.stdc.stdio;
import std.algorithm;
import std.array;
import std.string;
import mpi;
import mpi.util;

//#include <stdio.h>
//#include <stdlib.h>
//#include <mpi.h>
//#include <math.h>

/************************************************************
This is a simple broadcast program in MPI
************************************************************/

int main(string[] args)
{
    int argc = cast(int)args.length;
    auto argv = args.toArgv;

    int i,myid, numprocs;
    int source,count;
    int[4] buffer;
    MPI_Status status;
    MPI_Request request;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    source=0;
    count=4;
    if(myid == source)
    {
      for(i=0;i<count;i++)
        buffer[i]=i;
    }
    MPI_Bcast(buffer.ptr, count, MPI_INT, source, MPI_COMM_WORLD);

    for(i=0;i<count;i++)
      printf("%d ",buffer[i]);

    printf("\n");
    return MPI_Finalize();
}

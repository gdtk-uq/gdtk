///Shamelessly copied from:
///http://dvbmonkey.wordpress.com/2009/03/02/an-open-mpi-master-servant-example/
///If it destroys your system, blame that dude ^^


import std.stdio;
import std.array;
import std.algorithm;
import std.string;
import mpi;
import mpi.util;
//#include <unistd.h>

int main(string[] args)
{   //Converts D-style args to C-style.
    //pretty decent, since MPI will pass this along to the other processes
    int argc = cast(int)args.length;
    auto argv = args.toArgv(); 

    int numprocs, rank, namelen;
    char[MPI_MAX_PROCESSOR_NAME] processor_name;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name.ptr, &namelen);

    if (rank == 0)
    {
      printf("[%02d/%02d %s]: I am the master\n", rank, numprocs, processor_name.ptr );
      // Tell the servants to do something
    } else
    {
      printf("[%02d/%02d %s]: I am a servant\n", rank, numprocs, processor_name.ptr );
      // Wait for something to do
    }

    return MPI_Finalize();
}

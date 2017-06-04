///Shamelessly copied from:
///http://dvbmonkey.wordpress.com/2009/03/02/an-open-mpi-master-servant-example/
///If it destroys your system, blame that dude ^^
///No warranty, even the implied warrenty of usefulness or safety.

import std.stdio;
import core.stdc.stdio;
import core.stdc.string;
import std.array;
import std.algorithm;
import std.string;
import mpi;
import mpi.util;

int main(string[] args)
{
    int argc = cast(int)args.length;
    auto argv = args.toArgv();

    char[128] idstr;
    char[128] buff;
    char[256] stuff;


    int numprocs, rank, namelen, i;
    char[MPI_MAX_PROCESSOR_NAME] processor_name;

    MPI_Status stat;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name.ptr, &namelen);

    // Based on example from https://wiki.inf.ed.ac.uk/pub/ANC/ComputationalResources/slides.pdf
    if (rank == 0)
    {
        // This is the rank-0 copy of the process
        printf("Master Processor %d Reporting!\n", rank);
        printf("We have %d processors\n", numprocs);
        // Send each process a "Hello ... " string
        for(i = 1; i < numprocs; i++)
        {
            //Do we have sprintf??
            sprintf(buff.ptr, "Hello %d... ", i);
            MPI_Send(buff.ptr, 128, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        }
        // Go into a blocking-receive for each servant process
        for(i = 1; i < numprocs; i++)
        {
            MPI_Recv(buff.ptr, 128, MPI_CHAR, i, 0, MPI_COMM_WORLD, &stat);
            puts(buff.ptr);
        }
    }
    else
    {
        // Go into a blocking-receive waiting
        MPI_Recv(buff.ptr, 128, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &stat);
        // Append our identity onto the received string
        sprintf(idstr.ptr, "Processor %d Reporting!", rank);

        //Now here there should be some magic to concatenate the two
        //strings.  however, printf, sprintf and such do \0, and D
        //has shit support for that.
        //So instead of both, you only get one.
        //strcat is easy way to get around that.
        strcat(buff.ptr, idstr.ptr);

        // Send the string back to the rank-0 process
        MPI_Send(buff.ptr, 128, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
   }

   return MPI_Finalize();
}

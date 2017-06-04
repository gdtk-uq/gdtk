//Courtesy of Wikipedia and "My Lawyer!!"

import std.stdio;
import core.stdc.stdio;
import std.array;
import std.algorithm;
import std.string;
import core.stdc.string;
import mpi;
import mpi.util;

int main(string[] args)
{
    int argc = cast(int)args.length;
    auto argv = args.toArgv();

    string idstr;
    char[1024] buff;
    buff[] = 0;  // strings are inited with UTF invalid value by default
    size_t lastBuffPos;  // track last position

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
        writefln("Master Processor %d Reporting!", rank);
        writefln("We have %d processors", numprocs);

        // Send each process a "Hello ... " string
        for(i = 1; i < numprocs; i++)
        {
            auto str = format("Hello %s! ", i);
            buff[lastBuffPos .. lastBuffPos + str.length] = str[];
//            lastBuffPos += str.length;
            MPI_Send(buff.ptr, buff.length, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        }

        // Go into a blocking-receive for each servant process
        for(i = 1; i < numprocs; i++)
        {
            MPI_Recv(buff.ptr, buff.length, MPI_CHAR, i, 0, MPI_COMM_WORLD, &stat);
            writefln("%s: %s", rank, buff);
        }
    }
    else
    {
        // Go into a blocking-receive waiting
        MPI_Recv(buff.ptr, buff.length, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &stat);
        lastBuffPos = strlen(buff.ptr);  // calculate incoming string length

        // Append our identity onto the received string
        auto str = format("Processor %d reporting for duty!", rank);
        buff[lastBuffPos .. lastBuffPos + str.length] = str[];
        lastBuffPos += str.length;

        // Send the string back to the rank-0 process
        MPI_Send(buff.ptr, buff.length, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
   }

   return MPI_Finalize();
}
/** @Istrystfrf About that code: you can replace that stack-allocation with
 * heap-allocation, it doesn't matter. What matters is that you can't send data
 * of one length and receive data of a different length. If you use concatenation
 * then .length will change, and this is what you pass to MPI
 */

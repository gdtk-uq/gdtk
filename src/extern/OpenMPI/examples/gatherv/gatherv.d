///Found on http://geco.mines.edu/workshop/class2/examples/mpi/index.html

import std.stdio;
import core.stdc.stdlib;
import std.algorithm;
import std.array;
import std.string;
import core.memory;
import mpi;
import mpi.util;
/*
! This program shows how to use MPI_Gatherv.  Each processor sends a
! different amount of data to the root processor.  We use MPI_Gather
! first to tell the root how much data is going to be sent.
*/
/* globals */
int numnodes, myid, mpi_err;
enum mpi_root = 0;
/* end of globals */

void init_it(int* argc, char*** argv) {
    mpi_err = MPI_Init(argc, argv);
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &numnodes);
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}

void main(string[] args)
{
    int argc = cast(int)args.length;
    auto argv = args.toArgv();

    int* will_use;
    int* displacements, counts, allray;
    int size;

    init_it(&argc,&argv);
    auto mysize = myid + 1;
    auto myray = cast(int*)malloc(mysize * int.sizeof);
    myray[0 .. mysize] = mysize;

    /* counts and displacement arrays are only required on the root */
    if(myid == mpi_root){
        counts = cast(int*)malloc(numnodes * int.sizeof);
        displacements = cast(int*)malloc(numnodes * int.sizeof);
    }
    /* we gather the counts to the root */
    mpi_err = MPI_Gather(
            cast(void*)myray, 1, MPI_INT, 
            cast(void*)counts, 1, MPI_INT, 
            mpi_root, MPI_COMM_WORLD);
    /* calculate displacements and the size of the recv array */
    if(myid == mpi_root){
        displacements[0]=0;
        foreach(i; 1 .. numnodes){
            displacements[i] = counts[i-1] + displacements[i-1];
        }
        size=0;
        foreach(i; 0 .. numnodes)
            size = size + counts[i];
        allray = cast(int*)malloc(size * int.sizeof);
    }
    /* different amounts of data from each processor  */
    /* is gathered to the root */
    mpi_err = MPI_Gatherv(
            myray, mysize, MPI_INT, 
            allray, counts, displacements, MPI_INT,  
            mpi_root, MPI_COMM_WORLD);
                    
    if(myid == mpi_root){
        writeln(allray[0 .. size]);
    }
    mpi_err = MPI_Finalize();
}

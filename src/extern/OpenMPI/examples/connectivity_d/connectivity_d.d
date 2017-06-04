/*
 * Copyright (c) 2007      Sun Microsystems, Inc.  All rights reserved.
 */

/*
 * Test the connectivity between all processes.
 */

import std.stdio;
import std.algorithm;
import std.array;
//#include <errno.h>
//#include <stdio.h>
//#include <stdlib.h>
import core.stdc.stdlib;
//#include <string.h>
import std.string;
//#include <netdb.h>
//#include <unistd.h>
//#include <mpi.h>
import mpi;
import mpi.util;

int main(string[] args)
{
    int argc = cast(int)args.length;
    auto argv = args.toArgv();

    MPI_Status  status;
    bool        verbose = false;
    int         rank;
    int         np;        /* number of processes in job */
    int         peer;
    int         i;
    int         j;
    int         length;
    char[MPI_MAX_PROCESSOR_NAME + 1]        name;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    /*
     * If we cannot get the name for whatever reason, just
     * set it to unknown. */
    if (MPI_SUCCESS != MPI_Get_processor_name(name.ptr, &length)) {
        name[0.."unknown".length] = "unknown";
    }

    if (argc>1 && (cmp(args[1], "-v")==0))
        verbose = true;

    for (i=0; i<np; i++) {
        if (rank==i) {
            /* rank i sends to and receives from each higher rank */
            for(j=i+1; j<np; j++) {
                if (verbose)
                    printf("checking connection between rank %i on %s and rank %-4i\n", i, name.ptr, j);
                MPI_Send(&rank, 1, MPI_INT, j, rank, MPI_COMM_WORLD);
                MPI_Recv(&peer, 1, MPI_INT, j, j, MPI_COMM_WORLD, &status);
            }
        } else if (rank>i) {
            /* receive from and reply to rank i */
            MPI_Recv(&peer, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Send(&rank, 1, MPI_INT, i, rank, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0)
        printf("Connectivity test on %d processes PASSED.\n", np);

    MPI_Finalize();
    return 0;
}

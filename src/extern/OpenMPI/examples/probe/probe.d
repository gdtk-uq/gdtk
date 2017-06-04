/* Shows how to use probe & get_count to find the size of an */
/* incomming message                                         */
///Found on http://geco.mines.edu/workshop/class2/examples/mpi/index.html

import core.stdc.stdio;
import std.stdio;
import core.stdc.stdlib;
import std.string;
import std.algorithm;
import std.array;
import mpi;
import mpi.util;

void main(string[] args)
{
    int argc = cast(int)args.length;
    auto argv = args.toArgv();

    int myid, numprocs;
    MPI_Status status;
    int mytag,ierr,icount,j;
    int *i;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    printf(" Hello from D process: %d  Numprocs is %d\n",myid,numprocs);
    
    mytag=123;
    if(myid == 0) {
        j=200;
        icount=1;
        ierr=MPI_Send(&j,icount,MPI_INT,1,mytag,MPI_COMM_WORLD);
    }
    if(myid == 1){
        ierr=MPI_Probe(0,mytag,MPI_COMM_WORLD,&status);
        ierr=MPI_Get_count(&status,MPI_INT,&icount);
        i=cast(int*)malloc(icount * int.sizeof);
        printf("getting %d\n",icount);
        ierr = MPI_Recv(i,icount,MPI_INT,0,mytag,MPI_COMM_WORLD,&status);
        printf("i= ");
        for(j=0;j<icount;j++)
            printf("%d ",i[j]);
        printf("\n");
    }
    MPI_Finalize();
}

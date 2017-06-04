///Found on http://geco.mines.edu/workshop/class2/examples/mpi/index.html
//TODO newCommunicator doesn't work, check it out...

/****
! This program is designed to show how to set up a new communicator.
! We set up a communicator that includes all but one of the processors,
! The last processor is not part of the new communcator, TIMS_COMM_WORLD
! We use the routine MPI_Group_rank to find the rank within the new
! connunicator.  For the last processor the rank is MPI_UNDEFINED because
! it is not part of the communicator.  For this processor we call get_input
! The processors in TIMS_COMM_WORLD pass a token between themselves in the
! subroutine pass_token.  The remaining processor gets input, i, from the terminal
! and passes it to processor 1 of MPI_COMM_WORLD.  If i > 100 the program stops
*****/
import core.stdc.stdio;
import core.stdc.stdlib;
import std.algorithm;
import std.array;
import std.string;
import core.memory;
import mpi;
import mpi.util;

/* global variables */
int numnodes, myid, mpi_err;
enum mpi_root = 0;
/* end of global variables  */


void init_it(int* argc, char*** argv)
{
    mpi_err = MPI_Init(argc,argv);
    mpi_err = MPI_Comm_size( MPI_COMM_WORLD, &numnodes );
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}

int main(string[] args)
{
    int argc = cast(int)args.length;
    auto argv = args.toArgv();

/* poe a.out -procs 3 -rmpool 1 */
    int* will_use;
    MPI_Comm TIMS_COMM_WORLD;
    MPI_Group new_group, old_group;
    int num_used, used_id;

    init_it(&argc, &argv);
    printf("hello from %d\n", myid);

/* num_used is the # of processors that are part of the new communicator */
/* for this case hardwire to not include 1 processor */
    num_used = numnodes-1;
    if(numnodes > num_used)
    {
/* get our old group from MPI_COMM_WORLD */
        mpi_err = MPI_Comm_group(MPI_COMM_WORLD, &old_group);
/* create a new group from the old group */
/* that will contain a subset of the  processors */
        will_use = cast(int*)malloc(num_used * int.sizeof);
        foreach(ijk; 0 .. num_used)
        {
            will_use[ijk] = ijk;
        }
        mpi_err =  MPI_Group_incl(old_group, num_used, will_use, &new_group);
/* create the new communicator */
        mpi_err =  MPI_Comm_create(MPI_COMM_WORLD, new_group, &TIMS_COMM_WORLD);
/* test to see if I am part of new_group. */
        mpi_err =  MPI_Group_rank(new_group, &used_id);
        if(used_id == MPI_UNDEFINED)
        {
/* if not part of the new group do keyboard i/o. */
            get_input();
            mpi_err =  MPI_Barrier(MPI_COMM_WORLD);
            mpi_err =  MPI_Finalize();
            return mpi_err;
        }
    }
    else
    {
        mpi_err = MPI_Comm_dup(MPI_COMM_WORLD, &TIMS_COMM_WORLD );
    }
    mpi_err = MPI_Comm_rank(TIMS_COMM_WORLD, &myid );
    mpi_err = MPI_Comm_size(TIMS_COMM_WORLD, &numnodes );
    pass_token(TIMS_COMM_WORLD);
    printf("back\n");
    mpi_err = MPI_Barrier(MPI_COMM_WORLD);
    return MPI_Finalize();
}



void get_input()
{
   int to, my_tag, i;
   FILE *IN;
   to=0;
   my_tag=0;
   i=0;
   IN=fopen("ex10.in","r");
   assert(IN);
   while(i < 100)
   {
        printf("waiting for input:\n");
        fscanf(IN, "%i", &i);
        printf("i = %d\n",i);
        mpi_err = MPI_Send(cast(void*)&i, 1, MPI_INT, to, my_tag, MPI_COMM_WORLD);
   }
}

void pass_token(MPI_Comm THE_COMM_WORLD)
{
    int my_tag, j, i, to, from, ierr;
    MPI_Status status;
    int flag;
    my_tag = 1234;
    j = 0;
    flag = 1;
    to = myid+1;
    if(to == numnodes)
        to = 0;
    from = myid-1;
    if(from < 0)
        from = numnodes-1;
    i = 0;
    while(i < 100)
    {
        if(myid == j)
        {
            ierr = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG,
                    MPI_COMM_WORLD, &flag,&status);
            if(flag)
            {
                ierr = MPI_Recv(cast(void*)&i, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
                        MPI_COMM_WORLD, &status);
                printf("got i %d\n",i);
            }
            ierr = MPI_Send(cast(void*)&i, 1, MPI_INT, to, my_tag, THE_COMM_WORLD);
            j = -1;
        }
        else
        {
            ierr = MPI_Recv(cast(void*)&i, 1, MPI_INT, from, my_tag, THE_COMM_WORLD, &status);
            j = myid;
        }
    }
    if(myid < numnodes-1)
        ierr = MPI_Send(cast(void*)&i, 1, MPI_INT, to, my_tag, THE_COMM_WORLD);
}

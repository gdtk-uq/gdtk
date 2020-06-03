---
title: "Install & Run On HPC"
date: 2018-02-26T13:33:39+10:00
draft: false
---

# Install and Run on High-Performance Computers

This page describes how to install and run Eilmer on
some of the HPC cluster systems our research group
typically has access to.
If you have access to another cluster, you still might
find the information here useful.
You might find your local cluster setup is similar
to one of those listed here.
Likely, you should be able to adapt the instructions with
minimal hassle.

This section assumes you know how to set up your environment
for running eilmer on the cluster machine. What we give
you here are specific steps required to build and run
the MPI version.

{{< hint info >}}
A note on D compiler versions

As D is a relatively new language, and particularly so
in HPC environments, we find it best to bring-your-own
D compiler. We recommend having a very recent version
of the D compiler installed in your own account and
commit to updating the compiler on a regular basis.
For optimised builds use the LLVM compiler.
Install notes for the LLVM D compiler (LDC) are
available [here]({{< relref installing-ldc >}}).
{{< /hint >}}

## Goliath: EAIT Faculty Cluster

*Hardware:* Dell servers, Intel CPUs, 24-cores per node

*Operating system:* CentOS Linux 7.4

*Queue system:* SLURM

### Compiling
As this is a modern RedHat-flavoured system, you will need
to load the openmpi module to both build and run `e4mpi`.

To compile:

    module load mpi/openmpi-x86_64
    cd dgd/src/eilmer
    make FLAVOUR=fast WITH_MPI=1 install

### Running

Where possible, try to use a full node at a time on Goliath.
That means trying to set up your simulations to use 24 cores
per node.
Let's look at an example SLURM submit script where you
have 24 MPI tasks to run on 1 node.
I've named this file `run-on-goliath.slurm`.

    #!/bin/bash
    #SBATCH --job-name=my-MPI-job
    #SBATCH --nodes=1
    #SBATCH --ntasks=24

    module load mpi/openmpi-x86_64
    mpirun e4mpi --job=myJob --run  > LOGFILE

We submit that job using `sbatch`:

    sbatch run-on-goliath.slurm

Next, we'll look at running an oversubscribed job.
An oversubscribed job is when we setup *more* MPI tasks
than we have cores.
This might be the case if we have many blocks (eg. 120) and we
wish to run the 24 cores of a single node.


    #!/bin/bash
    #SBATCH --job-name=my-oversubscribed-MPI-job
    #SBATCH --nodes=1
    #SBATCH --ntasks=120
    #SBATCH --overcommit

    module load mpi/openmpi-x86_64
    mpirun e4mpi --job=myJob --run  > LOGFILE

Note the main changes are to correctly specify the number of
MPI tasks (here we have 120) and to add the directive `--overcommit`.


## Tinaroo: UQ RCC Cluster

*Hardware:* SGI severs, Intel CPUs, 24 cores-per-node, InfiniBand interconnect

*Operating system:* Rocks using CentOS Linux 6.9

*Queue system:* PBS

### Compiling
Tinaroo has an extensive module system with various options for
compilers and MPI libraries.
These have been compiled specifically on tinaroo to take
advantage of its InfiniBand interconnect.
Select an openmpi environment for compiling.
The `openmpi2_ib/2.0.2` has been tested and works nicely with the code.
To compile Eilmer4 with MPI enabled on tinaroo, do:

    module purge
    module load gnu
    module load openmpi2_ib/2.0.2
    cd dgd/src/install-scripts
    ./install-transient-solvers.sh

### Running
Where possible, try to use a full node at once on tinaroo.
The following is a tinaroo submit script that will request
24 MPI tasks. We have set up our job such that it is split
across 24 MPI tasks.

    #!/bin/bash
    #PBS -S /bin/bash
    #PBS -N my-MPI-job
    #PBS -A UQ-EAIT-MechMining
    #PBS -l select=1:ncpus=24:mpiprocs=24:mem=10g
    #PBS -l walltime=24:00:00

    module purge
    module load gnu
    module load openmpi2_ib/2.0.2

    cd $PBS_O_WORKDIR
    mpirun e4mpi --job=dbl-cone --run > LOGFILE

The `select` value is used to select the number of nodes.
The remaining values related to values per node.
So we've requested 24 CPUs on a node and 10 GB of memory
on a node.
If you set `select=3`, then the total resource request
would be for 72 CPUs and 30 GB of memory.
We use `qsub` to submit our job script:

    qsub run-on-tinaroo.qsub

## Gadi: NCI Cluster

*Hardware:* Fujitsu servers, Intel 'Cascade Lake' processors: 48 cores-per-node (two 24-core CPUs), 192 Gigabytes of RAM per node, HDR InfiniBand interconnect

*Operating system:* Linux, CentOS 8

*Queue system:* PBSPro

### An example setup in `.bash_profile`
Here, I include what I have added to the end of my `.bash_profile` file on gadi.
It sets my environment up to compile and run Eilmer.
Note also that I configure access to a locally installed version of the ldc2 compiler.

    export DGD=${HOME}/dgdinst
    append_path PATH ${DGD}/bin
    append_path PATH ${HOME}/opt/ldc2/bin
    append_path PYTHONPATH ${DGD}/lib
     
    export DGD_LUA_PATH=$DGD/lib/?.lua
    export DGD_LUA_CPATH=$DGD/lib/?.so
    
    module load openmpi/4.0.2

### Compiling

On gadi, you will need the ldc2 compiler installed locally.
In terms of modules, you only need to load the openmpi module.
I have had success with the 4.0.2 version.

    module load openmpi/4.0.2
    
As described for other systems, use the `install-transient-solvers.sh` script to
get an optimised build of the distributed-memory (MPI) transient solver.

    cd dgd/src/install-scripts
    ./install-transient-solvers.sh

### Running

As is common on large cluster computers, you will need to request the entire node CPU resources if your job spans multiple nodes. On gadi, that means CPU request numbers are in multiples of 48.
Here is a submit script where I've set up 192 MPI tasks for my job.

    #!/bin/bash
    #PBS -N my-MPI-job
    #PBS -P dc7
    #PBS -l walltime=00:30:00
    #PBS -l ncpus=192
    #PBS -l mem=200GB
    #PBS -l wd
    #PBS -l storage=scratch/dc7

    mpirun e4mpi --job=myJob --run > LOGFILE

Jobs on gadi are submitted using `qsub`.

Take note about the `storage` directive that appears in the submission script.
It is strongly encouraged on gadi to set up your jobs in the `scratch` area associated with your
project.
On gadi, you need to request explicit access to the `scratch` area in your submission
script so that filesystem is available to your job.
The `storage` directive used to request that access.

In the next example, I have set my job up to run in an oversubscribed mode.
Here I have 80 MPI tasks but I'd like to use only 48 CPUs.

    #!/bin/bash
    #PBS -N my-oversubscribed-MPI-job
    #PBS -P dc7
    #PBS -l walltime=00:30:00
    #PBS -l ncpus=48
    #PBS -l mem=50GB
    #PBS -l wd
    #PBS -l storage=scratch/dc7

    mpirun --oversubscribe -n 80 e4mpi --job=myJob --run > LOGFILE-oversubscribed



## NSCC Cluster in Singapore

*Hardware:* Fujitsu servers, Intel E5-2690v3 CPUs, 24 cores and 128 GB memory per node

*Operating system:* Red Hat Enterprise Linux Server release 6.9

*Queue system:* PBS

### Compiling

If git is not installed on this cluster, you may either install
it yourself, or just synch the source code across
from elsewhere. The second option is probably easiest. Use 'rsync'
to copy the source code from your own machine onto the cluster.

To load the modules required for compiling and running, and set up $PATH,
add the following to the `.bashrc` file:

    module purge
    module load openmpi/gcc493
    module load readline

    export DGD=$HOME/dgdinst
    export DGD_REPO=$HOME/dgd
    export PATH=$PATH:$DGD/bin
    export DGD_LUA_PATH=$DGD/lib/?.lua
    export DGD_LUA_CPATH=$DGD/lib/?.so

At the time of writing (16th May 2018), the openmpi module for gcc
installed on NSCC is for gcc4.9.3. You might need to check for newer versions
if it becomes obsolete in the future. Another thing to note is that once
the openmpi/gcc493 module is loaded, it automatically loads gcc/4.9.3 and
an extra bunch of C libraries.

To install the MPI version:

    cd dgd/src/eilmer
    make FLAVOUR=fast WITH_MPI=1 install

### Running

Where possible, try to use a full nodes at a time on NSCC. That means trying to set up your
simulations to use 24 cores per node. Here is a job script that uses 3 nodes or 72 CPUs.

    #!/bin/bash
    #PBS -q normal
    #PBS -N MyJob
    #PBS -l select=3:ncpus=24:mem=72G:mpiprocs=24:ompthreads=1
    #PBS -l walltime=23:59:59

    cd $PBS_O_WORKDIR
    mpirun e4mpi --job=MyJob --run > LogFile


In this example, LogFile contains the log of Eilmer execution,
 which can provide useful info for status checking and debugging.


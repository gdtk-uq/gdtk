/**
 * Machinery for solving solving linear systems in parallel
 *
 * Author: Nick Gibbons
 * Version: 2021-08-03: Prototyping
 */

module efieldexchange;

import std.stdio;
import std.math;
import std.algorithm;
import std.conv;

import fvinterface;
import geom;
import efieldbc;

version(mpi_parallel){
import mpi;

struct ExchangeData{
    int sendtag, recvtag;
    int nsend, nrecv;
    int sendrank, recvrank;
    double[] sendbuffer;
    double[] recvbuffer;
    int[] sendidxs, requestidxs;
    MPI_Request sendrequest;
    MPI_Status recvstatus;
}

class Exchanger {
/*
    Framework for doing sparse matrix operations over MPI. We take info
    from the block boundaries to figure out which processes have matrix
    entries we need. 

    Notes:
     - 
    @author: Nick Gibbons
*/
    double[] external_cell_buffer;

    this(const int rank, const int worldsize, const MPISharedField[] bcs) {
        this.rank = rank;
        this.worldsize = worldsize;

        // At boundaries, there is a pair of send/receive operations with the adjoining block
        foreach(bc; bcs){
            ExchangeData exch = ExchangeData();

            int nexch = to!int(bc.other_cell_ids.length);
            exch.sendtag = make_mpi_tag(rank, bc.other_blk_rank, worldsize);
            exch.recvtag = make_mpi_tag(bc.other_blk_rank, rank, worldsize);
            exch.nsend = nexch;
            exch.nrecv = nexch;
            exch.sendrank = bc.other_blk_rank; // send to
            exch.recvrank = bc.other_blk_rank; // receive from
            exch.sendbuffer.length = nexch;    // For matrix data
            exch.recvbuffer.length = nexch;    // For matrix data
            exch.sendidxs.length = nexch;      // send these cells to the other process
            exch.requestidxs.length = nexch;   // used to setup the other block's sendidxs

            foreach(i, cidx; bc.other_cell_ids){
                exch.requestidxs[i] = cidx;
            }
            exchanges ~= exch;
        }

        // external_cell_buffer has all of the data from all of the boundaries laid out on it
        // The order is each boundary in order of creation, then within a boundary, the order
        // matches the order of faces in the fluid boundary condition object.
        size_t N = 0;
        foreach(bc; bcs) N += bc.other_cell_ids.length;
        external_cell_buffer.length = N;

        // We use the bc.other_cell_ids to figure out which cells from the other block we need
        // We then send these over MPI to that other block. Under the hood, other_cell_ids
        // is assembled from the full_face_copy object in the boundary condition.
        foreach(exch; exchanges){
            MPI_Isend(exch.requestidxs.ptr, exch.nsend, MPI_INT, exch.sendrank, exch.sendtag, MPI_COMM_WORLD, &(exch.sendrequest));
        }

        // The other block will want some of our cells, sendidxs will contain these.
        foreach(exch; exchanges){
            MPI_Recv(exch.sendidxs.ptr, exch.nrecv, MPI_INT, exch.recvrank, exch.recvtag, MPI_COMM_WORLD, &(exch.recvstatus));
        }
    }

    void update_buffers(double[] arr){
    /*
        When performing a sparse matrix vector product in parallel, we need access to
        elements of the vector that are stored in another process. This routine copies
        those elements into external_cell_buffer, and sends some of our own over MPI
        to the other processes.

        Inputs:
         - arr : The vector being matrix multiplied

        Notes:
         - we could save some lines here by making each exch.recvbuffer a slice of
           external cell buffer.

        @author: Nick Gibbons
    */
        foreach(exch; exchanges){
            foreach(i, idx; exch.sendidxs){
                exch.sendbuffer[i] = arr[idx];
            }
            MPI_Isend(exch.sendbuffer.ptr, exch.nsend, MPI_DOUBLE, exch.sendrank, exch.sendtag, MPI_COMM_WORLD, &(exch.sendrequest));
        }

        size_t j=0;
        foreach(exch; exchanges){
            MPI_Recv(exch.recvbuffer.ptr, exch.nrecv, MPI_DOUBLE, exch.recvrank, exch.recvtag, MPI_COMM_WORLD, &(exch.recvstatus));
            foreach(recv; exch.recvbuffer){
                external_cell_buffer[j] = recv;
                j+=1;
            }
        }
    }

private:
    int rank, worldsize; 
    ExchangeData[] exchanges;

    @nogc
    int make_mpi_tag(int send_id, int recv_id, int worldsize){
    /*
        When doing a lagre number of many to many MPI operations, we need to
        tag each transfer with a unique identifier so they can be kept track of.

        Inputs:
         - send_id   : The MPI rank of the the sending process.
         - recv_id   : The MPI rank of the recevigin process.
         - worldsize : The total number of processes in the MPI job.

        Ouputs:
         - tag : A unique tag that can be used to arrange MPI transfers.

        Example usage for a 0->1 send:
         - Process 0 calls:
           tag = make_mpi_tag(0, 1, 2);
           or
           tag = make_mpi_rag(myrank, otherrank, 2);

         - Process 1 calls:
           tag = make_mpi_tag(0, 1, 2);
           or
           tag = make_mpi_rag(otherrank, myrank, 2);
    */
        return send_id + worldsize*recv_id;
    }
}
}

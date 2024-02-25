/**
 * Machinery for solving the large sparse matrix associated with a Parabolic Field Problem.
 *
 * Author: Nick Gibbons
 * Version: 2021-05-24: Prototyping
 */

module efieldgmres;

import std.stdio;
import std.math;
import std.algorithm;

import fvinterface;
import geom;
import efieldbc;
import efieldexchange;
version(mpi_parallel){
    import mpi;
}

class GMResFieldSolver {
    this() {}

    version(mpi_parallel){
        this(Exchanger exchanger) {
            this.exchanger = exchanger;
        }
    }

    void givens_rotation_cs(int i, int j, int n, double[] h, ref double c, ref double s){
        double a = h[j*n + j];
        double b = h[i*n + j];
        double r = sqrt(a*a + b*b);
        c = a/r;
        s = -b/r;
    }
    
    // Note that the caller has responsibility to allocate xd, unlike in the python code.
    void qr_lsq_backward_substitution(int n, int Rdim1, double[] R, double[] b, ref double[] xd){
        xd[] = 0.0;
        for (int k=n-1; k>=0; k--){
            double sum = 0.0;
            for (int i=n-1; i>k; i--){
                sum += R[k*Rdim1+i]*xd[i];
            }
    
            if (fabs(R[k*Rdim1+k])<1e-16){
                xd[k] = 0.0; // Numpy seems to do this to catch singularities. Why though?
            } else {
                xd[k] = (b[k] - sum)/R[k*Rdim1+k];
            }
        }
    }
    
    double dot_product(double[] a, double[] b, int n){
        double sum = 0.0;
        foreach(i; 0 .. n) sum += a[i]*b[i];
        version(mpi_parallel){
            MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        return sum;
    }
    
    double vector_norm(double[] a, int n){
        double sum = 0.0;
        foreach(i; 0 .. n) sum += a[i]*a[i];
        version(mpi_parallel){
            MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        return sqrt(sum);
    }
    
    void banded_matrix_vector_product(int n, int nb, double[] A, int[] Ai, double[] x, ref double[] b){
        b[] = 0.0;

        version(mpi_parallel){ exchanger.update_buffers(x);}

        foreach(i; 0 .. n){
            foreach(j; 0 .. nb){
                int jj = Ai[i*nb+j];
                if (jj<0) continue;

                double xjj;
                version(mpi_parallel){
                    if (jj>=n){
                        xjj = exchanger.external_cell_buffer[jj-n];
                    } else {
                        xjj = x[jj];
                    }
                } else {
                    xjj = x[jj];
                }
                b[i] += A[i*nb+j]*xjj;
            }
        }
    }
    
    void transpose_and_matrix_multiply(int n, int m, double[] A, double[] b, double[] x){
        /*
            Compute x=(A.T).dot(b), where A is a matrix and b is a vector.
    
            Inputs:
             n : Number of rows in matrix A
             m : Number of columns in matrix A
             A : The matrix to be multiplied (n,m)
             b : The vector to be multiplied (n) <- Note the dimension is n!
    
            Outputs:
             x : The result vector (m)
    
            Notes:
             - Like the other routines in this module, we assume the caller is responsible for allocating x 
             - This routine is never called along the distributed axis of the main matrix, so it doesn't
               need to be MPI aware, I think.
        */
    
        // This is the outer loop, iterating over the long axis, the columns of A, typically matrix_size
        for(int i=0; i<m; i++){ 
            // This is the inner loop, iterating over the short axis, the rows of A, typically k+1 in size
            x[i] = 0.0;
            for(int j=0; j<n; j++){
                x[i] += A[j*m + i]*b[j];
            }
        }
    }
    
    void solve(int matrix_size, int nbands, double[] A, int[] Ai, double[] b, double[] x0,
                      double[] xf, int nmax_iter, bool verbose, double tol=1e-12){
        /*
            Generalised Minimal Residual Method for solving linear systems representated by a banded matrix.
    
            Based on: 
             - https://en.wikipedia.org/wiki/Generalized_minimal_residual_method
             - https://en.wikipedia.org/wiki/QR_decomposition
    
            Inputs:
             matrix_size : Number of rows in the matrix being solved for
             nbands      : Number of entries in a band, typically 5 for this application
             A           : Matrix to be solved for  (matrix_size,nbands)
             Ai          : Connectivity matrix specifying where each entry in Ai would be in the full matrix
             b           : RHS matrix (matrix_size)
             x0          : initial guess to begin building Krylov vectors (matrix_size)
             nmax_iter   : Maximum number of iterations to allocate space for
             tol         : Stop when residual falls below this number
    
            Outputs:
             xf          : The final answer vector (matrix_size)
    
            @author: Nick Gibbons
        */
        double[] r;
        r.length = matrix_size;
        banded_matrix_vector_product(matrix_size, nbands, A, Ai, x0, r);
        foreach(i; 0 .. matrix_size) r[i] = b[i] - r[i];
        double rnorm = vector_norm(r, matrix_size);
    
        double[] q,y,xold,xnew,xdiff,h,c,s,QT,R,B,Y;
    
        // Memory allocation 👀
        q.length   = nmax_iter*matrix_size;
        y.length   = matrix_size;
        xold.length= matrix_size;
        xnew.length= matrix_size;
        xdiff.length= matrix_size;

        h.length   = (nmax_iter+1)*nmax_iter;
        c.length   = nmax_iter;
        s.length   = nmax_iter;
        QT.length  = (nmax_iter+1)*(nmax_iter+1);
        R.length   = (nmax_iter+1)*nmax_iter;
        B.length   = nmax_iter+1;
        Y.length   = nmax_iter+1;
    
        QT[] = 0.0;
        foreach(i; 0 .. nmax_iter+1) QT[i*(nmax_iter+1)+i] = 1.0;
        R[] = 0.0;
        xold[] = x0[];
    
        // First basis vector is created from r
        foreach(i; 0 .. matrix_size) q[i] = r[i]/rnorm;
    
        int niters = nmax_iter;
        int nprint = nmax_iter/40;
        double residual = 1e99;
        bool success = true;
        int k;
    
        for(k=0; k<niters; k++){
            // Perform Arnoldi Iteration to generate basis vectors
            double[] qk = q[k*matrix_size .. (k+1)*matrix_size];
            banded_matrix_vector_product(matrix_size, nbands, A, Ai, qk, y);

            //writeln("qk");
            //writeln(q[k*matrix_size .. (k+1)*matrix_size]);
    
            for(int j; j<k+1; j++){
                double[] qj = q[j*matrix_size .. (j+1)*matrix_size];
                double hjk = dot_product(y, qj, matrix_size);
                h[j*nmax_iter + k] = hjk;
                foreach(p; 0 .. matrix_size) y[p] -= hjk*qj[p];
            }
    
            double hkp1k = vector_norm(y, matrix_size);
            h[(k+1)*nmax_iter + k] = hkp1k;
            if ((hkp1k != 0.0) && (k != nmax_iter-1)){
                foreach(p; 0 .. matrix_size) q[(k+1)*matrix_size + p] = y[p]/hkp1k;
            }
    
            // Now we need to grow the QR matrices that are used for the least squares problem
            int m = k+2;
            int n = k+1;
            // Start by copying the new column of h into R, the upper triangular matrix
            for(int i=0; i<m; i++) R[i*nmax_iter + k] = h[i*nmax_iter + k];
    
            // Now replay the old givens rotations across the new column, excluding the last element
            for(int j=0; j<k; j++){
                int i = j+1;
                double Rjk = R[j*nmax_iter + k];
                double Rik = R[i*nmax_iter + k];
                R[j*nmax_iter + k] = c[j]*Rjk + -s[j]*Rik;
                R[i*nmax_iter + k] = s[j]*Rjk +  c[j]*Rik;
                double Qjn = QT[j*(nmax_iter+1) + n];
                double Qin = QT[i*(nmax_iter+1) + n];
                QT[i*(nmax_iter+1) + n] = s[j]*Qjn + c[j]*Qin;
            }
    
            // Next, compute the new Givens rotation from the last element, and store it for replaying
            //writeln("h");
            //foreach(i; 0 .. k+2) writeln(h[i*nmax_iter .. (i*nmax_iter+k+1)]);
            int j = k;
            int i = k+1;
            givens_rotation_cs(i, j, nmax_iter, R, c[k], s[k]);
            double Rjk = R[j*nmax_iter + k];
            double Rik = R[i*nmax_iter + k];
            R[j*nmax_iter + k] = c[j]*Rjk + -s[j]*Rik;
            R[i*nmax_iter + k] = s[j]*Rjk +  c[j]*Rik;
            for(int p=0; p<m; p++){
                double Qjp = QT[j*(nmax_iter+1) + p];
                double Qip = QT[i*(nmax_iter+1) + p];
                QT[j*(nmax_iter+1) + p] = c[j]*Qjp + -s[j]*Qip;
                QT[i*(nmax_iter+1) + p] = s[j]*Qjp +  c[j]*Qip;
            }
    
            //writeln("QT");
            //foreach(p; 0 .. m) writeln(QT[p*(nmax_iter+1) .. p*(nmax_iter+1)+m]);
            //writeln("R");
            //foreach(p; 0 .. n) writeln(R[p*(nmax_iter) .. p*(nmax_iter)+n]);
    
            // With everything set up, solve the least squared problem and get a new answer vector xnew
            for(int p=0; p<n; p++) B[p] = rnorm*QT[p*(nmax_iter+1) + 0];
            qr_lsq_backward_substitution(n, nmax_iter, R, B, Y);
            transpose_and_matrix_multiply(n, matrix_size, q, Y, xnew);
            for(int p=0; p<matrix_size; p++) xnew[p] += x0[p];

            for(int p=0; p<matrix_size; p++) xdiff[p] = (xnew[p] - xold[p]);
            residual = vector_norm(xdiff, matrix_size);
            xold[] = xnew[];
    
            //if ((k%nprint==0) && verbose) write(".");
            //writefln("iter: %d residual %e hkp1k %e", k, residual, hkp1k);
            if (residual<tol) break;
        }
        if (residual>=tol) success=false;
        if (verbose) writefln("    Solve Complete: status=%s  iters=%d/%d  residual=%e/%e", success, k, nmax_iter, residual, tol);
        if (success==false) throw new Error("BGMRes failed to converge!");

        xf[] = xnew[];
        return;
    }

private:
    version(mpi_parallel){
        Exchanger exchanger;
    }
}

void test_bgmres(){
    int nmax_iter = 5;
    int matrix_size = 5;
    int nbands = 3;

    double[] A = [ 0., -2.,  1.,
                   1.,  4.,  1.,
                   1., -5.,  1.,
                   1.,  5.,  1.,
                   1., -3.,  0.];
          
    int[] Ai = [-1, 0, 1,
                 0, 1, 2,
                 1, 2, 3,
                 2, 3, 4,
                 3, 4,-1];

    auto gmres = new GMResFieldSolver();
    double[] xtarget = [-1.0, 2.0, -1.2, 3.0, 0.8];
    double[] b;
    b.length = 5;
    gmres.banded_matrix_vector_product(5, 3, A, Ai, xtarget, b);
	double[] x;
    double[] x0;
    x.length = 5;
    x0.length = 5;

    x0[] = 0.0;

    gmres.solve(matrix_size, nbands, A, Ai, b, x0, x, nmax_iter, true);
    double error=0.0;
    foreach(p; 0 .. xtarget.length) error = (xtarget[p] - x[p])*(xtarget[p] - x[p]);
    error = sqrt(error);
    writeln("x:");
    writeln(x);
    writeln("xtarget:");
    writeln(xtarget);
    writeln("error: ", error);
}

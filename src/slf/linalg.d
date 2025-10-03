// Linear algebra for solving block sparse tridiagonal systems
// for slf, @author: NNG

module slf.linalg;

import std.conv;
import std.format;
import std.stdio;
import std.string;

import slf.io;

import nm.bbla;
import nm.number;
import ntypes.complex;

void multiply_tridigonal_block_matrix(ref const Parameters pm, Matrix!(double)[3][] LHS, Matrix!(double)[] U, Matrix!(double)[] R){
    size_t neq = pm.neq;
    size_t N   = pm.N;


    R[0] = LHS[0][1].dot(U[0]) + LHS[0][2].dot(U[0+1]);
    foreach(cell; 1 .. N-1){
        R[cell] = LHS[cell][0].dot(U[cell-1]) + LHS[cell][1].dot(U[cell]) + LHS[cell][2].dot(U[cell+1]);
    }
    R[N-1] = LHS[N-1][0].dot(U[N-2]) + LHS[N-1][1].dot(U[N-1]);
    return;

}

void solve_tridigonal_block_matrix(ref const Parameters pm, ref TridiagonalSolveWorkspace tdws, Matrix!(double)[3][] LHS, Matrix!(double)[] U, Matrix!(double)[] R){
/*
    Solve a sparse tridiagonal block matrix problem. We assume the data is packed as follows:
    [[##, B0, C0],
     [A1, B1, C1],
         ...
     [AN-1, BN-1, ##]]

*/

    size_t neq = pm.neq;
    size_t N   = pm.N;
    auto Ai = tdws.Ai;
    auto Bi = tdws.Bi;
    auto Ci = tdws.Ci;
    auto Ui = tdws.Ui;
    auto Ri = tdws.Ri;

    // First eliminate the A blocks using block row multiply and adds
    foreach(step; 0 .. N-1){
        Bi._data[] = LHS[step][1]._data[];
        Ci._data[] = LHS[step][2]._data[];
        Ri._data[] = R[step]._data[];

        auto perm = decomp!double(Bi);
        solve!double(Bi, Ri, perm);
        solve!double(Bi, Ci, perm);
        R[step]._data[] = Ri._data[];
        LHS[step][2]._data[] = Ci._data[];
        LHS[step][1].eye();

        Ai._data[] = LHS[step+1][0]._data[];
        LHS[step+1][0].zeros();
        LHS[step+1][1] -= Ai.dot(Ci);
        R[step+1]      -= Ai.dot(Ri);
    }

    // Now do a backward substituion on the remaining blocks to solve for U
    size_t end = N-1;
    Bi._data[] = LHS[end][1]._data[];
    Ri._data[] = R[end]._data[];

    auto perm = decomp!double(Bi);
    solve!double(Bi, Ri, perm);
    U[end]._data[] = Ri._data[];

    // Danger Looping down an unsigned integer is a bad idea...
    for (size_t step=N-2; step>=0; step--){
        Bi._data[] = LHS[step][1]._data[];
        Ci._data[] = LHS[step][2]._data[];
        Ri._data[] = R[step]._data[];
        Ui._data[] = U[step+1]._data[];

        Ri -= Ci.dot(Ui);

        perm = decomp!double(Bi);
        solve!double(Bi, Ri, perm);
        U[step]._data[] = Ri._data[];
        if (step==0) break;
    }

    return;
}





struct TridiagonalSolveWorkspace{
    Matrix!double Ai;
    Matrix!double Bi;
    Matrix!double Ci;
    Matrix!double Ui;
    Matrix!double Ri;

    this(size_t neq){
        Ai = new Matrix!double(neq, neq);
        Bi = new Matrix!double(neq, neq);
        Ci = new Matrix!double(neq, neq);
        Ui  = new Matrix!double(neq);
        Ri  = new Matrix!double(neq);
    }
}


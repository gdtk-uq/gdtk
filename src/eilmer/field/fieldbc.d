/**
 * fieldbcs.d
 * Machinery for applying boundary conditions to electromagnetic fields
 *
 * Author: Nick Gibbons
 * Version: 2021-05-24: Prototyping
 */

module fieldbc;

import std.stdio;
import std.math;
import std.format;
import std.conv;
import std.json;

import fvcell;
import fvinterface;
import geom;
import json_helper;
import nm.number;
import fieldconductivity;
import bc.boundary_condition;
import bc.ghost_cell_effect.full_face_copy;

interface FieldBC {
    void opCall(const double sign, const FVInterface face, const FVCell cell, ref double Akk, ref double bk, ref double Ako, ref int Aio);
    double compute_current(const double sign, const FVInterface face, const FVCell cell);
}

class ZeroNormalGradient : FieldBC {
    this() {}

    final void opCall(const double sign, const FVInterface face, const FVCell cell, ref double Akk, ref double bk, ref double Ako, ref int Aio){
        bk  = 0.0;
        Akk = 0.0;
        Ako = 0.0;
        Aio = -1;
    }
    final double compute_current(const double sign, const FVInterface face, const FVCell cell){
        return 0.0;
    }
}

class FixedField : FieldBC {
    this(double value) {
        this.value = value;
    }

    final void opCall(const double sign, const FVInterface face, const FVCell cell, ref double Akk, ref double bk, ref double Ako, ref int Aio){
        double S = face.length.re;
        double d = distance_between(face.pos, cell.pos[0]);
        double sigma = face.fs.gas.sigma.re;
        double phi = value;

        bk  = -1.0*phi*S/d*sigma;
        Akk = -1.0*S/d*sigma;
        Ako = 0.0;
        Aio = -1;
    }
    final double compute_current(const double sign, const FVInterface face, const FVCell cell){
        double S = face.length.re;
        double d = distance_between(face.pos, cell.pos[0]);
        double phigrad = (value - cell.electric_potential)/d; // This implicitly points out of the domain.
        double sigma = face.fs.gas.sigma.re;
        double I = sigma*phigrad*S; // convert from current density vector j to current I
        return I;
    }
private:
    double value;
}


class FixedField_Test : FieldBC {
    this() {
    }

    final void opCall(const double sign, const FVInterface face, const FVCell cell, ref double Akk, ref double bk, ref double Ako, ref int Aio){
        double S = face.length.re;
        double d = distance_between(face.pos, cell.pos[0]);
        double sigma = face.fs.gas.sigma.re;
        double phi = test_field(face.pos.x.re, face.pos.y.re);

        bk  = -1.0*phi*S/d*sigma;
        Akk = -1.0*S/d*sigma;
        Ako = 0.0;
        Aio = -1;
    }
    final double compute_current(const double sign, const FVInterface face, const FVCell cell){
        double S = face.length.re;
        double d = distance_between(face.pos, cell.pos[0]);
        double phi = test_field(face.pos.x.re, face.pos.y.re);
        double phigrad = (phi - cell.electric_potential)/d; // This implicitly points out of the domain.
        double sigma = face.fs.gas.sigma.re;
        double I = sigma*phigrad*S;
        return I;
    }
private:

    double test_field(double x, double y){
        return exp(x)*sin(y);
    }
}

class FixedGradient_Test : FieldBC {
    this() {
    }

    final void opCall(const double sign, const FVInterface face, const FVCell cell, ref double Akk, ref double bk, ref double Ako, ref int Aio){
        double S = face.length.re;
        double d = distance_between(face.pos, cell.pos[0]);
        double sigma = face.fs.gas.sigma.re;
        Vector3 phigrad = test_field_gradient(face.pos.x.re, face.pos.y.re);

        number phigrad_dot_n = phigrad.dot(face.n);
        bk  = -sign*phigrad_dot_n.re*S*sigma;
        Akk = 0.0;
        Ako = 0.0;
        Aio = -1;
    }
    final double compute_current(const double sign, const FVInterface face, const FVCell cell){
        double S = face.length.re;
        Vector3 phigrad = test_field_gradient(face.pos.x.re, face.pos.y.re);
        number phigrad_dot_n = sign*phigrad.dot(face.n); // TODO: Should this be negative sign?
        double sigma = face.fs.gas.sigma.re;
        double I = sigma*phigrad_dot_n.re*S;
        return I;
    }
private:

    Vector3 test_field_gradient(double x, double y){
        Vector3 phigrad = Vector3(exp(x)*sin(y), exp(x)*cos(y), 0.0);
        return phigrad;
    }
}

class SharedField : FieldBC {
    this(const BoundaryCondition bc, const int[] block_offsets) {
        GhostCellFullFaceCopy gc;
        foreach(action; bc.preReconAction){
            gc = cast(GhostCellFullFaceCopy) action;
            if (gc !is null){
                break;
            }
        }
        if (gc is null){
            throw new Error("Boundary is missing full_face_copy. Possibly check all outer boundaires are assigned.");
        }

        // We keep our own copies of data related to the shared boundary
        other_blk_id = gc.other_blk.id;
        other_block_offset = block_offsets[other_blk_id];

        // Since the arrays inside the GhostCellFullFaceCopy are not in the same order as the boundary faces
        // we have to do some work to organise our own mapping array, of boundary faces to shared cell ids.
        other_cell_ids.length = bc.faces.length;
        foreach(i, f; bc.faces){
            foreach(j, c; gc.ghost_cells){
                if ((f.right_cell==c) || (f.left_cell==c)) {
                    other_cell_ids[i] = to!int(gc.mapped_cell_ids[j]);
                    break;
                }
            }
        }
        return;
    }

    final void opCall(const double sign, const FVInterface face, const FVCell cell, ref double Akk, ref double bk, ref double Ako, ref int Aio){
        Vector3 ghost_cell_position;

        if (cell==face.left_cell){
            ghost_cell_position = face.right_cell.pos[0];
        } else {
            ghost_cell_position = face.left_cell.pos[0];
        }

        double S = face.length.re;
        double sigma = face.fs.gas.sigma.re;
        double d = distance_between(ghost_cell_position, cell.pos[0]);

        Aio = other_cell_ids[face.i_bndry] + other_block_offset; // CHECKME
        Ako =  1.0*S/d*sigma;
        Akk = -1.0*S/d*sigma;
        bk = 0.0;
    }
    final double compute_current(const double sign, const FVInterface face, const FVCell cell){
        return 0.0;
    }
private:
    int other_blk_id;
    int other_block_offset;
    int[] other_cell_ids;
}

version(mpi_parallel){
class MPISharedField : FieldBC {
    static int nExtraCells=0;
    static MPISharedField[] instances;
    int other_blk_rank;
    int other_blk_face;
    int other_blk_id;
    int[] other_cell_ids;

    this(const BoundaryCondition bc, const int ncells) {

        GhostCellFullFaceCopy gc;
        foreach(action; bc.preReconAction){
            gc = cast(GhostCellFullFaceCopy) action;
            if (gc !is null){
                break;
            }
        }
        if (gc is null){
            throw new Error("Boundary is missing full_face_copy. Possibly check all outer boundaries are assigned.");
        }

        // In MPI mode we assume that each block has its own process. We need to know:
        // 1.) Where in the extra cell array each othercell goes, for setting Ai in the matrix
        // 2.) We need to know where in the other block that other cell goes, so that the sending
        //     process can send it to us.
        // This extra cell array will be order of faces, but also in order of the SharedBCs being initialised.
        other_blk_rank = gc.other_blk_rank;
        other_blk_id = gc.other_blk.id;
        other_blk_face = gc.other_face;

        // Since the arrays inside the GhostCellFullFaceCopy are not in the same order as the boundary faces
        // we have to do some work to organise our own mapping array, of boundary faces to shared cell ids.
        other_cell_ids.length = bc.faces.length;
        external_cell_idxs.length = bc.faces.length;
        my_offset = nExtraCells + ncells;

        foreach(i, f; bc.faces){
            foreach(j, c; gc.ghost_cells){
                if ((f.right_cell==c) || (f.left_cell==c)) {
                    other_cell_ids[i] = to!int(gc.mapped_cell_ids[j]);
                    external_cell_idxs[i] = to!int(i) + my_offset;
                    break;
                }
            }
        }
        nExtraCells += bc.faces.length;
        instances ~= this;
        return;
    }

    final void opCall(const double sign, const FVInterface face, const FVCell cell, ref double Akk, ref double bk, ref double Ako, ref int Aio){
        Vector3 ghost_cell_position;

        if (cell==face.left_cell){
            ghost_cell_position = face.right_cell.pos[0];
        } else {
            ghost_cell_position = face.left_cell.pos[0];
        }

        double S = face.length.re;
        double sigma = face.fs.gas.sigma.re;
        double d = distance_between(ghost_cell_position, cell.pos[0]);

        Aio = external_cell_idxs[face.i_bndry];
        Ako =  1.0*S/d*sigma;
        Akk = -1.0*S/d*sigma;
        bk = 0.0;
    }
    final double compute_current(const double sign, const FVInterface face, const FVCell cell){
        return 0.0;
    }
private:
    int[] external_cell_idxs;
    int my_offset;
} // end class MPISharedField
} // end version(mpi_parallel)

FieldBC create_field_bc(JSONValue field_bc_json, const BoundaryCondition bc, const int[] block_offsets, string conductivity_model_name, int ncells){
/*
    Create a field_bc object that will be used later for setting matrix entries near boundaries.
    Currently the data specifying each bc is storied in bc.field_bc as a JSON table.

*/
    string name = getJSONstring(field_bc_json, "name", "not_found");
    FieldBC field_bc;

    switch (name) {
    case "ZeroNormalGradient":
        field_bc = new ZeroNormalGradient();
        break;
    case "FixedField":
        double value = getJSONdouble(field_bc_json, "value", 0.0);
        field_bc = new FixedField(value);
        break;
    case "FixedGradient_Test":
        field_bc = new FixedGradient_Test();
        break;
    case "FixedField_Test":
        field_bc = new FixedField_Test();
        break;
    case "unspecified":
        version(mpi_parallel){
            field_bc = new MPISharedField(bc, ncells);
        } else {
            field_bc = new SharedField(bc, block_offsets);
        }
        break;
    default:
        string errMsg = format("The FieldBC '%s' is not available.", name);
        throw new Error(errMsg);
    }
    return field_bc;
}

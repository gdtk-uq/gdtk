// lsqinterp.d
// Least-squares interpolation/reconstruction of flow field.
//

import std.math;
import std.stdio;
import nm.rsla;
import geom;
import gas;
import fvcore;
import globalconfig;
import flowstate;
import fvinterface;
import fvcell;

class LsqInterpolator {

private:
    LocalConfig myConfig;

public:
    this(LocalConfig myConfig) 
    {
	this.myConfig = myConfig;
    }

    @nogc double clip_to_limits(double q, double A, double B)
    // Returns q if q is between the values A and B, else
    // it returns the closer limit of the range [A,B].
    {
	double lower_limit = fmin(A, B);
	double upper_limit = fmax(A, B);
	return fmin(upper_limit, fmax(lower_limit, q));
    } // end clip_to_limits()

    void interp_both(ref FVInterface IFace, size_t gtl, ref FlowState Lft, ref FlowState Rght)
    {
	auto gmodel = myConfig.gmodel;
	auto nsp = gmodel.n_species;
	auto nmodes = gmodel.n_modes;
	// Low-order reconstruction just copies data from adjacent FV_Cell.
	// Even for high-order reconstruction, we depend upon this copy for
	// the viscous-transport and diffusion coefficients.
	Lft.copy_values_from(IFace.left_cells[0].fs);
	Rght.copy_values_from(IFace.right_cells[0].fs);
	if ( myConfig.interpolation_order > 1 ) {
	    // High-order reconstruction for some properties.
	    //
	    // Always reconstruct in the interface-local frame of reference.
	    foreach (cell; IFace.left_cells) {
		cell.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
	    }
	    foreach (cell; IFace.right_cells) {
		cell.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
	    }
	    // Prepare the normal matrix for each cloud and invert it.
	    // Since we are working in the interface-local frame, having the
	    // local x-direction aligned with the unit normal for the interface,
	    // we always expect good distribition of points in that direction.
	    // Depending on how aligned the points are, there may be insufficient
	    // variation in the local y- or z-positions to make a sensible 3D matrix.
	    // We will form the sums and look at their relative sizes.
	    Vector3 dr;
	    double[6][3] xTxL; // normal matrix Left, augmented to give 6 entries per row
	    double[3] rhsL, gradientsL;
	    size_t nL = IFace.left_cells.length;
	    double xx = 0.0; double xy = 0.0; double xz = 0.0;
	    double yy = 0.0; double yz = 0.0;
	    double zz = 0.0;
	    foreach (i; 1 .. nL) {
		dr = IFace.left_cells[i].pos[gtl] - IFace.left_cells[0].pos[gtl];
		dr.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
		xx += dr.x*dr.x; xy += dr.x*dr.y; xz += dr.x*dr.z;
		yy += dr.y*dr.y; yz += dr.y*dr.z; zz += dr.z*dr.z;
	    }
	    assert (fabs(xx) > 1.0e-12, "left_cells xx essentially zero");
	    string directionsL = "x";
	    if (yy/xx > 0.01) directionsL ~= "y";
	    if (zz/xx > 0.01) directionsL ~= "z";
	    switch (directionsL) {
	    case "xyz":
		xTxL[0][0] = xx; xTxL[0][1] = xy; xTxL[0][2] = xz;
		xTxL[1][0] = xy; xTxL[1][1] = yy; xTxL[1][2] = yz;
		xTxL[2][0] = xz; xTxL[2][1] = yz; xTxL[2][2] = zz;
		break;
	    case "xy":
		xTxL[0][0] =  xx; xTxL[0][1] =  xy; xTxL[0][2] = 0.0;
		xTxL[1][0] =  xy; xTxL[1][1] =  yy; xTxL[1][2] = 0.0;
		xTxL[2][0] = 0.0; xTxL[2][1] = 0.0; xTxL[2][2] = 1.0;
		break;
	    case "xz":
		xTxL[0][0] =  xx; xTxL[0][1] = 0.0; xTxL[0][2] =  xz;
		xTxL[1][0] = 0.0; xTxL[1][1] = 1.0; xTxL[1][2] = 0.0;
		xTxL[2][0] =  xz; xTxL[2][1] = 0.0; xTxL[2][2] =  zz;
		break;
	    case "x":
	    default:
		xTxL[0][0] =  xx; xTxL[0][1] = 0.0; xTxL[0][2] = 0.0;
		xTxL[1][0] = 0.0; xTxL[1][1] = 1.0; xTxL[1][2] = 0.0;
		xTxL[2][0] = 0.0; xTxL[2][1] = 0.0; xTxL[2][2] = 1.0;
	    }
	    xTxL[0][3] = 1.0; xTxL[0][4] = 0.0; xTxL[0][5] = 0.0;
	    xTxL[1][3] = 0.0; xTxL[1][4] = 1.0; xTxL[1][5] = 0.0;
	    xTxL[2][3] = 0.0; xTxL[2][4] = 0.0; xTxL[2][5] = 1.0;
	    if (0 != computeInverse!3(xTxL)) {
		// Assume that the rows are linearly dependent 
		// because the sample points are colinear.
		// Proceed by working as a single-dimensional interpolation.
		directionsL = "x";
		xTxL[0][0] = 1.0; xTxL[0][1] = 0.0; xTxL[0][2] = 0.0;
		xTxL[1][0] = 0.0; xTxL[1][1] = 1.0; xTxL[1][2] = 0.0;
		xTxL[2][0] = 0.0; xTxL[2][1] = 0.0; xTxL[2][2] = 1.0;
		xTxL[0][3] = 1.0/xx; xTxL[0][4] = 0.0; xTxL[0][5] = 0.0;
		xTxL[1][3] = 0.0; xTxL[1][4] = 1.0; xTxL[1][5] = 0.0;
		xTxL[2][3] = 0.0; xTxL[2][4] = 0.0; xTxL[2][5] = 1.0;
	    }
	    //
	    double[6][3] xTxR; // normal matrix Right, augmented to give 6 entries per row
	    double[3] rhsR, gradientsR;
	    size_t nR = IFace.right_cells.length;
	    xx = 0.0; xy = 0.0; xz = 0.0;
	    yy = 0.0; yz = 0.0;
	    zz = 0.0;
	    foreach (i; 1 .. nR) {
		dr = IFace.right_cells[i].pos[gtl] - IFace.right_cells[0].pos[gtl];
		dr.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
		xx += dr.x*dr.x; xy += dr.x*dr.y; xz += dr.x*dr.z;
		yy += dr.y*dr.y; yz += dr.y*dr.z; zz += dr.z*dr.z;
	    }
	    assert (fabs(xx) > 1.0e-12, "right_cells xx essentially zero");
	    string directionsR = "x";
	    if (yy/xx > 0.01) directionsR ~= "y";
	    if (zz/xx > 0.01) directionsR ~= "z";
	    switch (directionsR) {
	    case "xyz":
		xTxR[0][0] = xx; xTxR[0][1] = xy; xTxR[0][2] = xz;
		xTxR[1][0] = xy; xTxR[1][1] = yy; xTxR[1][2] = yz;
		xTxR[2][0] = xz; xTxR[2][1] = yz; xTxR[2][2] = zz;
		break;
	    case "xy":
		xTxR[0][0] =  xx; xTxR[0][1] =  xy; xTxR[0][2] = 0.0;
		xTxR[1][0] =  xy; xTxR[1][1] =  yy; xTxR[1][2] = 0.0;
		xTxR[2][0] = 0.0; xTxR[2][1] = 0.0; xTxR[2][2] = 1.0;
		break;
	    case "xz":
		xTxR[0][0] =  xx; xTxR[0][1] = 0.0; xTxR[0][2] =  xz;
		xTxR[1][0] = 0.0; xTxR[1][1] = 1.0; xTxR[1][2] = 0.0;
		xTxR[2][0] =  xz; xTxR[2][1] = 0.0; xTxR[2][2] =  zz;
		break;
	    case "x":
	    default:
		xTxR[0][0] =  xx; xTxR[0][1] = 0.0; xTxR[0][2] = 0.0;
		xTxR[1][0] = 0.0; xTxR[1][1] = 1.0; xTxR[1][2] = 0.0;
		xTxR[2][0] = 0.0; xTxR[2][1] = 0.0; xTxR[2][2] = 1.0;
	    }
	    xTxR[0][3] = 1.0; xTxR[0][4] = 0.0; xTxR[0][5] = 0.0;
	    xTxR[1][3] = 0.0; xTxR[1][4] = 1.0; xTxR[1][5] = 0.0;
	    xTxR[2][3] = 0.0; xTxR[2][4] = 0.0; xTxR[2][5] = 1.0;
	    if (0 != computeInverse!3(xTxR)) {
		// Assume that the rows are linearly dependent 
		// because the sample points are colinear.
		// Proceed by working as a single-dimensional interpolation.
		directionsR = "x";
		xTxR[0][0] = 1.0; xTxR[0][1] = 0.0; xTxR[0][2] = 0.0;
		xTxR[1][0] = 0.0; xTxR[1][1] = 1.0; xTxR[1][2] = 0.0;
		xTxR[2][0] = 0.0; xTxR[2][1] = 0.0; xTxR[2][2] = 1.0;
		xTxR[0][3] = 1.0/xx; xTxR[0][4] = 0.0; xTxR[0][5] = 0.0;
		xTxR[1][3] = 0.0; xTxR[1][4] = 1.0; xTxR[1][5] = 0.0;
		xTxR[2][3] = 0.0; xTxR[2][4] = 0.0; xTxR[2][5] = 1.0;
	    }
	    // x-velocity
	    string codeForReconstruction(string qname, string tname)
	    {
		string code = "
                foreach (j; 0 .. 3) { rhsL[j] = 0.0; }
                foreach (i; 1 .. nL) {
                    double dq = IFace.left_cells[i].fs."~qname~" - IFace.left_cells[0].fs."~qname~";
		    dr = IFace.left_cells[i].pos[gtl] - IFace.left_cells[0].pos[gtl];
		    dr.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    rhsL[0] += dr.x*dq; rhsL[1] += dr.y*dq; rhsL[2] += dr.z*dq;
	        }
	        solveWithInverse!3(xTxL, rhsL, gradientsL);
		switch (directionsL) {
		case `xyz`: break;
		case `xy`: gradientsL[2] = 0.0; break;
		case `xz`: gradientsL[1] = 0.0; break;
		case `x`: default: gradientsL[1] = 0.0; gradientsL[2] = 0.0;
                }
                foreach (j; 0 .. 3) { rhsR[j] = 0.0; }
                foreach (i; 1 .. nR) {
                    double dq = IFace.right_cells[i].fs."~qname~" - IFace.right_cells[0].fs."~qname~";
 		    dr = IFace.right_cells[i].pos[gtl] - IFace.right_cells[0].pos[gtl];
		    dr.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    rhsR[0] += dr.x*dq; rhsR[1] += dr.y*dq; rhsR[2] += dr.z*dq;
	        }
	        solveWithInverse!3(xTxR, rhsR, gradientsR);
		switch (directionsR) {
		case `xyz`: break;
		case `xy`: gradientsR[2] = 0.0; break;
		case `xz`: gradientsR[1] = 0.0; break;
		case `x`: default: gradientsR[1] = 0.0; gradientsR[2] = 0.0;
                }
		// TODO limiting of gradients.
                dr = IFace.pos - IFace.left_cells[0].pos[gtl];
		dr.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                Lft."~tname~" = IFace.left_cells[0].fs."~qname~";
                Lft."~tname~" += dr.x * gradientsL[0] + dr.y * gradientsL[1] + dr.z * gradientsL[2];
                Lft."~tname~" = clip_to_limits(Lft."~tname~", IFace.left_cells[0].fs."~qname~",
                                               IFace.right_cells[0].fs."~qname~");
                dr = IFace.pos - IFace.right_cells[0].pos[gtl];
		dr.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                Rght."~tname~" = IFace.right_cells[0].fs."~qname~";
                Rght."~tname~" += dr.x * gradientsR[0] + dr.y * gradientsR[1] + dr.z * gradientsR[2];
                Rght."~tname~" = clip_to_limits(Rght."~tname~", IFace.left_cells[0].fs."~qname~",
                                                IFace.right_cells[0].fs."~qname~");
                ";
		return code;
	    }
	    mixin(codeForReconstruction("vel.x", "vel.refx"));
	    mixin(codeForReconstruction("vel.y", "vel.refy"));
	    mixin(codeForReconstruction("vel.z", "vel.refz"));
	    if (myConfig.MHD) {
		mixin(codeForReconstruction("B.x", "B.refx"));
		mixin(codeForReconstruction("B.y", "B.refy"));
		mixin(codeForReconstruction("B.z", "B.refz"));
	    }
	    if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
		mixin(codeForReconstruction("tke", "tke"));
		mixin(codeForReconstruction("omega", "omega"));
	    }
	    if (nsp > 1) {
		// Multiple species.
		for ( size_t isp = 0; isp < nsp; ++isp ) {
		    mixin(codeForReconstruction("gas.massf[isp]", "gas.massf[isp]"));
		}
		try {
		    scale_mass_fractions(Lft.gas.massf); 
		} catch {
		    Lft.gas.massf[] = IFace.left_cells[0].fs.gas.massf[];
		}
		try {
		    scale_mass_fractions(Rght.gas.massf);
		} catch {
		    Rght.gas.massf[] = IFace.right_cells[0].fs.gas.massf[];
		}
	    } else {
		// Only one possible mass-fraction value for a single species.
		Lft.gas.massf[0] = 1.0;
		Rght.gas.massf[0] = 1.0;
	    }
	    // Interpolate on two of the thermodynamic quantities, 
	    // and fill in the rest based on an EOS call. 
	    // If an EOS call fails, fall back to just copying cell-centre data.
	    // This does presume that the cell-centre data is valid. 
	    final switch (myConfig.thermo_interpolator) {
	    case InterpolateOption.pt: 
		mixin(codeForReconstruction("gas.p", "gas.p"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstruction("gas.T[imode]", "gas.T[imode]"));
		}
		try {
		    gmodel.update_thermo_from_pT(Lft.gas);
		} catch {
		    Lft.copy_values_from(IFace.left_cells[0].fs);
		}
		try {
		    gmodel.update_thermo_from_pT(Rght.gas);
		} catch {
		    Rght.copy_values_from(IFace.right_cells[0].fs);
		}
		break;
	    case InterpolateOption.rhoe:
		mixin(codeForReconstruction("gas.rho", "gas.rho"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstruction("gas.e[imode]", "gas.e[imode]"));
		}
		try {
		    gmodel.update_thermo_from_rhoe(Lft.gas);
		} catch {
		    Lft.copy_values_from(IFace.left_cells[0].fs);
		}
		try {
		    gmodel.update_thermo_from_rhoe(Rght.gas);
		} catch {
		    Rght.copy_values_from(IFace.right_cells[0].fs);
		}
		break;
	    case InterpolateOption.rhop:
		mixin(codeForReconstruction("gas.rho", "gas.rho"));
		mixin(codeForReconstruction("gas.p", "gas.p"));
		try {
		    gmodel.update_thermo_from_rhop(Lft.gas);
		} catch {
		    Lft.copy_values_from(IFace.left_cells[0].fs);
		}
		try {
		    gmodel.update_thermo_from_rhop(Rght.gas);
		} catch {
		    Rght.copy_values_from(IFace.right_cells[0].fs);
		}
		break;
	    case InterpolateOption.rhot: 
		mixin(codeForReconstruction("gas.rho", "gas.rho"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstruction("gas.T[imode]", "gas.T[imode]"));
		}
		try {
		    gmodel.update_thermo_from_rhoT(Lft.gas);
		} catch {
		    Lft.copy_values_from(IFace.left_cells[0].fs);
		}
		try {
		    gmodel.update_thermo_from_rhoT(Rght.gas);
		} catch {
		    Rght.copy_values_from(IFace.right_cells[0].fs);
		}
		break;
	    } // end switch thermo_interpolator
	    // Finally, undo the transformation to local coordinates.
	    Lft.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
	    Rght.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
	    foreach (cell; IFace.left_cells) {
		cell.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
	    }
	    foreach (cell; IFace.right_cells) {
		cell.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
	    }
	} // end of high-order reconstruction
    } // end interp_both()

} // end class LsqInterpolator

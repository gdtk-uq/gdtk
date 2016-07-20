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


class LSQInterpWorkspace {
    // A place to hold the intermediate results for assembling
    // the normal equations, to allow reuse of the values.
public:
    double dxFace, dyFace, dzFace;
    double[] dx, dy, dz;
    double[6][3] xTx; // normal matrix, augmented to give 6 entries per row

    this()
    {
	// don't need to do anything
    }

    this(const LSQInterpWorkspace other)
    {
	dx = other.dx.dup(); dy = other.dy.dup(); dz = other.dz.dup();
	dxFace = other.dxFace; dyFace = other.dyFace; dzFace = other.dzFace;
	foreach (i; 0 .. 3) {
	    foreach (j; 0 .. 6) { xTx[i][j] = other.xTx[i][j]; }
	}
    }
} // end class LSQInterpWorkspace


class LsqInterpolator {

private:
    LocalConfig myConfig;
    LSQInterpWorkspace common_wsL, common_wsR;

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

    @nogc void min_mod_limit(ref double a, ref double b)
    // A rough slope limiter.
    {
	if (a * b < 0.0) {
	    a = 0.0; b = 0.0;
	    return;
	}
	a = copysign(fmin(fabs(a), fabs(b)), a);
	b = a;
    }

    @nogc void van_albada_limit(ref double a, ref double b)
    // A smooth slope limiter.
    {
	immutable double eps = 1.0e-12;
	double s = (a*b + fabs(a*b))/(a*a + b*b + eps);
	a *= s;
	b *= s;
    }

    void assemble_and_invert_normal_matrix(ref FVInterface IFace, size_t gtl,
					   ref BasicCell[] cell_cloud,
					   ref LSQInterpWorkspace ws)
    {
	// Since we are working in the interface-local frame, having the
	// local x-direction aligned with the unit normal for the interface,
	// we always expect good distribition of points in that direction.
	Vector3 dr = IFace.pos; dr -= cell_cloud[0].pos[gtl];
	dr.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
	ws.dxFace = dr.x; ws.dyFace = dr.y; ws.dzFace = dr.z;
	//
	size_t np = cell_cloud.length;
	if (ws.dx.length < np) { ws.dx.length = np; }
	if (ws.dy.length < np) { ws.dy.length = np; }
	if (ws.dz.length < np) { ws.dz.length = np; }
	foreach (i; 1 .. np) {
	    dr = cell_cloud[i].pos[gtl]; dr -= cell_cloud[0].pos[gtl];
	    dr.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
	    ws.dx[i] = dr.x; ws.dy[i] = dr.y; ws.dz[i] = dr.z;
	}	    
	// Prepare the normal matrix for the cloud and invert it.
	// Depending on how aligned the points are, there may be insufficient
	// variation in the local y- or z-positions to make a sensible 3D matrix.
	// We will form the sums and look at their relative sizes.
	double xx = 0.0; double xy = 0.0; double xz = 0.0;
	double yy = 0.0; double yz = 0.0;
	double zz = 0.0;
	foreach (i; 1 .. np) {
	    xx += ws.dx[i]*ws.dx[i]; xy += ws.dx[i]*ws.dy[i]; xz += ws.dx[i]*ws.dz[i];
	    yy += ws.dy[i]*ws.dy[i]; yz += ws.dy[i]*ws.dz[i]; zz += ws.dz[i]*ws.dz[i];
	}
	assert (fabs(xx) > 1.0e-12, "left_cells xx essentially zero");
	if (myConfig.dimensions == 3) {
	    ws.xTx[0][0] = xx; ws.xTx[0][1] = xy; ws.xTx[0][2] = xz;
	    ws.xTx[1][0] = xy; ws.xTx[1][1] = yy; ws.xTx[1][2] = yz;
	    ws.xTx[2][0] = xz; ws.xTx[2][1] = yz; ws.xTx[2][2] = zz;
	} else {
	    // dimensions == 2
	    ws.xTx[0][0] =  xx; ws.xTx[0][1] =  xy; ws.xTx[0][2] = 0.0;
	    ws.xTx[1][0] =  xy; ws.xTx[1][1] =  yy; ws.xTx[1][2] = 0.0;
	    ws.xTx[2][0] = 0.0; ws.xTx[2][1] = 0.0; ws.xTx[2][2] = 1.0;
	}
	ws.xTx[0][3] = 1.0; ws.xTx[0][4] = 0.0; ws.xTx[0][5] = 0.0;
	ws.xTx[1][3] = 0.0; ws.xTx[1][4] = 1.0; ws.xTx[1][5] = 0.0;
	ws.xTx[2][3] = 0.0; ws.xTx[2][4] = 0.0; ws.xTx[2][5] = 1.0;
	if (0 != computeInverse!(3,3)(ws.xTx)) {
	    // Assume that the rows are linearly dependent 
	    // because the sample points are colinear.
	    // Proceed by working as a single-dimensional interpolation.
	    ws.xTx[0][0] = 1.0; ws.xTx[0][1] = 0.0; ws.xTx[0][2] = 0.0;
	    ws.xTx[1][0] = 0.0; ws.xTx[1][1] = 1.0; ws.xTx[1][2] = 0.0;
	    ws.xTx[2][0] = 0.0; ws.xTx[2][1] = 0.0; ws.xTx[2][2] = 1.0;
	    ws.xTx[0][3] = 1.0/xx; ws.xTx[0][4] = 0.0; ws.xTx[0][5] = 0.0;
	    ws.xTx[1][3] = 0.0; ws.xTx[1][4] = 0.0; ws.xTx[1][5] = 0.0;
	    ws.xTx[2][3] = 0.0; ws.xTx[2][4] = 0.0; ws.xTx[2][5] = 0.0;
	}
    } // end assemble_and_invert_normal_matrix()
    
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
	if (myConfig.interpolation_order > 1) {
	    // High-order reconstruction for some properties.
	    //
	    LSQInterpWorkspace wsL;
	    LSQInterpWorkspace wsR;
	    if (myConfig.retain_least_squares_work_data) {
		// To have a faster calculation (by nearly a factor of 2),
		// retain the workspaces with precomputed inverses in the interfaces.
		wsL = IFace.wsL;
		wsR = IFace.wsR;
	    } else {
		// To reduce memory consumption, we use common workspaces.
		wsL = common_wsL;
		wsR = common_wsR;
		assemble_and_invert_normal_matrix(IFace, gtl, IFace.left_cells, wsL);
		assemble_and_invert_normal_matrix(IFace, gtl, IFace.right_cells, wsR);
	    }
	    //
	    // Always reconstruct in the interface-local frame of reference.
	    // note that clouds may share the interface neighbour cells, you should be sure not to
	    // transform the neighbour cells twice!
	    foreach (left_cell; IFace.left_cells) {
		if (left_cell.id != IFace.right_cells[0].id){
		    left_cell.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
		}
	    }
	    foreach (right_cell; IFace.right_cells) {
		if (right_cell.id != IFace.left_cells[0].id){
		    right_cell.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
		}
	    }
	    //
	    // Actual resonstruction phase.
	    double[3] rhsL, gradientsL;
	    double[3] rhsR, gradientsR;
	    // x-velocity
	    string codeForReconstruction(string qname, string tname)
	    {
		string code = "{
                double qL0 = IFace.left_cells[0].fs."~qname~";
                foreach (j; 0 .. 3) { rhsL[j] = 0.0; }
                foreach (i; 1 .. IFace.left_cells.length) {
                    double dq = IFace.left_cells[i].fs."~qname~" - qL0;
                    rhsL[0] += wsL.dx[i]*dq; rhsL[1] += wsL.dy[i]*dq; rhsL[2] += wsL.dz[i]*dq;
	        }
	        solveWithInverse!(3,3)(wsL.xTx, rhsL, gradientsL);
                double qR0 = IFace.right_cells[0].fs."~qname~";
                foreach (j; 0 .. 3) { rhsR[j] = 0.0; }
                foreach (i; 1 .. IFace.right_cells.length) {
                    double dq = IFace.right_cells[i].fs."~qname~" - qR0;
                    rhsR[0] += wsR.dx[i]*dq; rhsR[1] += wsR.dy[i]*dq; rhsR[2] += wsR.dz[i]*dq;
	        }
                solveWithInverse!(3,3)(wsR.xTx, rhsR, gradientsR);
                if (myConfig.apply_limiter) {
		    van_albada_limit(gradientsL[0], gradientsR[0]);
                    van_albada_limit(gradientsL[1], gradientsR[1]);
                }
                double qL = qL0 + wsL.dxFace * gradientsL[0] + wsL.dyFace * gradientsL[1];
                double qR = qR0 + wsR.dxFace * gradientsR[0] + wsR.dyFace * gradientsR[1];
                if (myConfig.dimensions == 3) {
                    if (myConfig.apply_limiter) { van_albada_limit(gradientsL[2], gradientsR[2]); }
                    qL += wsL.dzFace * gradientsL[2];
                    qR += wsR.dzFace * gradientsR[2];
                }
                if (myConfig.extrema_clipping) {
                    Lft."~tname~" = clip_to_limits(qL, qL0, qR0);
                    Rght."~tname~" = clip_to_limits(qR, qL0, qR0);
                } else {
                    Lft."~tname~" = qL;
                    Rght."~tname~" = qR;
                }
                }
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
		if (myConfig.divergence_cleaning) {
		    mixin(codeForReconstruction("psi", "psi"));
		}
	    }
	    if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
		mixin(codeForReconstruction("tke", "tke"));
		mixin(codeForReconstruction("omega", "omega"));
	    }
	    if (nsp > 1) {
		// Multiple species.
		foreach (isp; 0 .. nsp) {
		    mixin(codeForReconstruction("gas.massf[isp]", "gas.massf[isp]"));
		}
		try {
		    scale_mass_fractions(Lft.gas.massf); 
		} catch (Exception e) {
		    writeln(e.msg);
		    Lft.gas.massf[] = IFace.left_cells[0].fs.gas.massf[];
		}
		try {
		    scale_mass_fractions(Rght.gas.massf);
		} catch (Exception e) {
		    writeln(e.msg);
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
	    string codeForThermoUpdate(string funname)
	    {
		string code = "
		try {
		    gmodel.update_thermo_from_"~funname~"(Lft.gas);
		} catch (Exception e) {
		    writeln(e.msg);
		    Lft.copy_values_from(IFace.left_cells[0].fs);
		}
		try {
		    gmodel.update_thermo_from_"~funname~"(Rght.gas);
		} catch (Exception e) {
		    writeln(e.msg);
		    Rght.copy_values_from(IFace.right_cells[0].fs);
		}
                ";
		return code;
	    }
	    final switch (myConfig.thermo_interpolator) {
	    case InterpolateOption.pt: 
		mixin(codeForReconstruction("gas.p", "gas.p"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstruction("gas.T[imode]", "gas.T[imode]"));
		}
		mixin(codeForThermoUpdate("pT"));
		break;
	    case InterpolateOption.rhoe:
		mixin(codeForReconstruction("gas.rho", "gas.rho"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstruction("gas.e[imode]", "gas.e[imode]"));
		}
		mixin(codeForThermoUpdate("rhoe"));
		break;
	    case InterpolateOption.rhop:
		mixin(codeForReconstruction("gas.rho", "gas.rho"));
		mixin(codeForReconstruction("gas.p", "gas.p"));
		mixin(codeForThermoUpdate("rhop"));
		break;
	    case InterpolateOption.rhot: 
		mixin(codeForReconstruction("gas.rho", "gas.rho"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstruction("gas.T[imode]", "gas.T[imode]"));
		}
		mixin(codeForThermoUpdate("rhoT"));
		break;
	    } // end switch thermo_interpolator
	    // Finally, undo the transformation to local coordinates.
	    Lft.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
	    Rght.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
	    foreach (left_cell; IFace.left_cells) {
		if (left_cell.id != IFace.right_cells[0].id){
		    left_cell.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
		}
	    }
	    foreach (right_cell; IFace.right_cells) {
		if (right_cell.id != IFace.left_cells[0].id){
		    right_cell.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
		}
	    }
	} // end of high-order reconstruction
    } // end interp_both()
    
} // end class LsqInterpolator


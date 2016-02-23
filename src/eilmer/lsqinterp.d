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

// Symbolic names for the reconstruction directions.
enum Directions { x, xy, xz, xyz }

class LsqInterpolator {

private:
    LocalConfig myConfig;
    double dxFaceL, dyFaceL, dzFaceL;
    double dxFaceR, dyFaceR, dzFaceR;
    immutable size_t MaxNPoints = 10;
    double[MaxNPoints] dxL, dyL, dzL;
    double[MaxNPoints] dxR, dyR, dzR;

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
    // Try a bit of a rough limiter on the slopes.
    {
	if (a * b < 0.0) {
	    a = 0.0; b = 0.0;
	    return;
	}
	a = copysign(fmin(fabs(a), fabs(b)), a);
	b = a;
    }
    
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
	    // Since we are working in the interface-local frame, having the
	    // local x-direction aligned with the unit normal for the interface,
	    // we always expect good distribition of points in that direction.
	    Vector3 dr = IFace.pos; dr -= IFace.left_cells[0].pos[gtl];
	    dr.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
	    dxFaceL = dr.x; dyFaceL = dr.y; dzFaceL = dr.z;
	    dr = IFace.pos; dr -= IFace.right_cells[0].pos[gtl];
	    dr.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
	    dxFaceR = dr.x; dyFaceR = dr.y; dzFaceR = dr.z;
	    //
	    size_t nL = IFace.left_cells.length;
	    assert(nL <= MaxNPoints, "too many left_cells");
	    foreach (i; 1 .. nL) {
		dr = IFace.left_cells[i].pos[gtl]; dr -= IFace.left_cells[0].pos[gtl];
		dr.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
		dxL[i] = dr.x; dyL[i] = dr.y; dzL[i] = dr.z;
	    }	    
	    size_t nR = IFace.right_cells.length;
	    assert(nR < MaxNPoints, "too many right_cells");
	    foreach (i; 1 .. nR) {
		dr = IFace.right_cells[i].pos[gtl]; dr -= IFace.right_cells[0].pos[gtl];
		dr.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
		dxR[i] = dr.x; dyR[i] = dr.y; dzR[i] = dr.z;
	    }	    
	    // Prepare the normal matrix for each cloud and invert it.
	    // Depending on how aligned the points are, there may be insufficient
	    // variation in the local y- or z-positions to make a sensible 3D matrix.
	    // We will form the sums and look at their relative sizes.
	    double[6][3] xTxL; // normal matrix Left, augmented to give 6 entries per row
	    double[3] rhsL, gradientsL;
	    double xx = 0.0; double xy = 0.0; double xz = 0.0;
	    double yy = 0.0; double yz = 0.0;
	    double zz = 0.0;
	    foreach (i; 1 .. nL) {
		xx += dxL[i]*dxL[i]; xy += dxL[i]*dyL[i]; xz += dxL[i]*dzL[i];
		yy += dyL[i]*dyL[i]; yz += dyL[i]*dzL[i]; zz += dzL[i]*dzL[i];
	    }
	    assert (fabs(xx) > 1.0e-12, "left_cells xx essentially zero");
	    immutable double small = 0.01;
	    auto directionsL = Directions.xyz;
	    if (yy/xx > small && zz/xx > small) {
		directionsL = Directions.xyz;
	    } else if (yy/xx > small && zz/xx <= small) {
		directionsL = Directions.xy;
	    } else if (yy/xx <= small && zz/xx > small) {
		directionsL = Directions.xz;
	    } else {
		directionsL = Directions.x;
	    }
	    final switch (directionsL) {
	    case Directions.xyz:
		xTxL[0][0] = xx; xTxL[0][1] = xy; xTxL[0][2] = xz;
		xTxL[1][0] = xy; xTxL[1][1] = yy; xTxL[1][2] = yz;
		xTxL[2][0] = xz; xTxL[2][1] = yz; xTxL[2][2] = zz;
		break;
	    case Directions.xy:
		xTxL[0][0] =  xx; xTxL[0][1] =  xy; xTxL[0][2] = 0.0;
		xTxL[1][0] =  xy; xTxL[1][1] =  yy; xTxL[1][2] = 0.0;
		xTxL[2][0] = 0.0; xTxL[2][1] = 0.0; xTxL[2][2] = 1.0;
		break;
	    case Directions.xz:
		xTxL[0][0] =  xx; xTxL[0][1] = 0.0; xTxL[0][2] =  xz;
		xTxL[1][0] = 0.0; xTxL[1][1] = 1.0; xTxL[1][2] = 0.0;
		xTxL[2][0] =  xz; xTxL[2][1] = 0.0; xTxL[2][2] =  zz;
		break;
	    case Directions.x:
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
		directionsL = Directions.x;
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
	    xx = 0.0; xy = 0.0; xz = 0.0;
	    yy = 0.0; yz = 0.0;
	    zz = 0.0;
	    foreach (i; 1 .. nR) {
		xx += dxR[i]*dxR[i]; xy += dxR[i]*dyR[i]; xz += dxR[i]*dzR[i];
		yy += dyR[i]*dyR[i]; yz += dyR[i]*dzR[i]; zz += dzR[i]*dzR[i];
	    }
	    assert (fabs(xx) > 1.0e-12, "right_cells xx essentially zero");
	    auto directionsR = Directions.xyz;
	    if (yy/xx > small && zz/xx > small) {
		directionsR = Directions.xyz;
	    } else if (yy/xx > small && zz/xx <= small) {
		directionsR = Directions.xy;
	    } else if (yy/xx <= small && zz/xx > small) {
		directionsR = Directions.xz;
	    } else {
		directionsR = Directions.x;
	    }
	    final switch (directionsR) {
	    case Directions.xyz:
		xTxR[0][0] = xx; xTxR[0][1] = xy; xTxR[0][2] = xz;
		xTxR[1][0] = xy; xTxR[1][1] = yy; xTxR[1][2] = yz;
		xTxR[2][0] = xz; xTxR[2][1] = yz; xTxR[2][2] = zz;
		break;
	    case Directions.xy:
		xTxR[0][0] =  xx; xTxR[0][1] =  xy; xTxR[0][2] = 0.0;
		xTxR[1][0] =  xy; xTxR[1][1] =  yy; xTxR[1][2] = 0.0;
		xTxR[2][0] = 0.0; xTxR[2][1] = 0.0; xTxR[2][2] = 1.0;
		break;
	    case Directions.xz:
		xTxR[0][0] =  xx; xTxR[0][1] = 0.0; xTxR[0][2] =  xz;
		xTxR[1][0] = 0.0; xTxR[1][1] = 1.0; xTxR[1][2] = 0.0;
		xTxR[2][0] =  xz; xTxR[2][1] = 0.0; xTxR[2][2] =  zz;
		break;
	    case Directions.x:
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
		directionsR = Directions.x;
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
		string code = "{
                double qL0 = IFace.left_cells[0].fs."~qname~";
                foreach (j; 0 .. 3) { rhsL[j] = 0.0; }
                foreach (i; 1 .. nL) {
                    double dq = IFace.left_cells[i].fs."~qname~" - qL0;
                    rhsL[0] += dxL[i]*dq; rhsL[1] += dyL[i]*dq; rhsL[2] += dzL[i]*dq;
	        }
	        solveWithInverse!3(xTxL, rhsL, gradientsL);
		final switch (directionsL) {
		case Directions.xyz: break;
		case Directions.xy: gradientsL[2] = 0.0; break;
		case Directions.xz: gradientsL[1] = 0.0; break;
		case Directions.x: gradientsL[1] = 0.0; gradientsL[2] = 0.0;
                }
                double qR0 = IFace.right_cells[0].fs."~qname~";
                foreach (j; 0 .. 3) { rhsR[j] = 0.0; }
                foreach (i; 1 .. nR) {
                    double dq = IFace.right_cells[i].fs."~qname~" - qR0;
                    rhsR[0] += dxR[i]*dq; rhsR[1] += dyR[i]*dq; rhsR[2] += dzR[i]*dq;
	        }
                solveWithInverse!3(xTxR, rhsR, gradientsR);
		final switch (directionsR) {
		case Directions.xyz: break;
		case Directions.xy: gradientsR[2] = 0.0; break;
		case Directions.xz: gradientsR[1] = 0.0; break;
		case Directions.x: gradientsR[1] = 0.0; gradientsR[2] = 0.0;
                }
		min_mod_limit(gradientsL[0], gradientsR[0]);
                double qL = qL0 + dxFaceL * gradientsL[0] + 
                            dyFaceL * gradientsL[1] + dzFaceL * gradientsL[2];
                Lft."~tname~" = clip_to_limits(qL, qL0, qR0);
                double qR = qR0 + dxFaceR * gradientsR[0] + 
                            dyFaceR * gradientsR[1] + dzFaceR * gradientsR[2];
                Rght."~tname~" = clip_to_limits(qR, qL0, qR0);
                }";
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
	    final switch (myConfig.thermo_interpolator) {
	    case InterpolateOption.pt: 
		mixin(codeForReconstruction("gas.p", "gas.p"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstruction("gas.T[imode]", "gas.T[imode]"));
		}
		try {
		    gmodel.update_thermo_from_pT(Lft.gas);
		} catch (Exception e) {
		    writeln(e.msg);
		    Lft.copy_values_from(IFace.left_cells[0].fs);
		}
		try {
		    gmodel.update_thermo_from_pT(Rght.gas);
		} catch (Exception e) {
		    writeln(e.msg);
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
		} catch (Exception e) {
		    writeln(e.msg);
		    Lft.copy_values_from(IFace.left_cells[0].fs);
		}
		try {
		    gmodel.update_thermo_from_rhoe(Rght.gas);
		} catch (Exception e) {
		    writeln(e.msg);
		    Rght.copy_values_from(IFace.right_cells[0].fs);
		}
		break;
	    case InterpolateOption.rhop:
		mixin(codeForReconstruction("gas.rho", "gas.rho"));
		mixin(codeForReconstruction("gas.p", "gas.p"));
		try {
		    gmodel.update_thermo_from_rhop(Lft.gas);
		} catch (Exception e) {
		    writeln(e.msg);
		    Lft.copy_values_from(IFace.left_cells[0].fs);
		}
		try {
		    gmodel.update_thermo_from_rhop(Rght.gas);
		} catch (Exception e) {
		    writeln(e.msg);
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
		} catch (Exception e) {
		    writeln(e.msg);
		    Lft.copy_values_from(IFace.left_cells[0].fs);
		}
		try {
		    gmodel.update_thermo_from_rhoT(Rght.gas);
		} catch (Exception e) {
		    writeln(e.msg);
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

/** BDF1, BDF2 and variable BDF2
*develops matrix class
*solves with LU decomposition using Crout's method
*/


module nm.bdfLU;

import std.stdio;
import std.math;
import std.conv;
import std.algorithm;
import std.exception;
import nm.bbla;

//BDF1 or backwards euler formula
//applies newton-raphson method to solve
//requires current state vector, current rate vector, jacobian (dR/dS), and timestep
//returns new guess state vector for newton method

double[] BDF1(double dt, double[] state, double[] state_guess, double[] rate, double[][] jacobian) {
	//Solution takes the form Ax=b
	//[I-dt(dR/dS)](dS) = state + dt*rate - state_guess
	auto A = new Matrix!double(state.length,state.length);
	double[] x;
	x.length=state.length;
	double[] b;
	b.length=state.length;
	//populate A and b
	for (int i;i<state.length;i++) {
		for (int j;j<state.length;j++) {
			if (i==j) {
				A[i,j]=1-dt*jacobian[i][j];
			} else {
				A[i,j]=-dt*jacobian[i][j];
			}
		}
		b[i]=state[i] + dt*rate[i] - state_guess[i];
	}
	//solve system with bbla
	int[] pivot;
	pivot.length=state.length;
	LUDecomp!double(A,pivot);
	LUSolve!double(A,pivot,b,x);
	
	double[] update_guess1=state_guess.dup;
	update_guess1[] += x[];
	
	return update_guess1;
}

double[] BDF2(double dt, double[] prev_state, double[] state, double[] state_guess, double[] rate, double[][] jacobian) {
	//Solution takes the form Ax=b
	//[I-(2/3)*dt(dR/dS)](dS) = (4/3)*state + (2/3)*dt*rate - (1/3)*prev_state - state_guess
	auto A = new Matrix!double(state.length,state.length);
	double[] x;
	x.length=state.length;
	double[] b;
	b.length=state.length;
	//populate A and b
	for (int i;i<state.length;i++) {
		for (int j;j<state.length;j++) {
			if (i==j) {
				A[i,j]=1-(2.0/3.0)*dt*jacobian[i][j];
			} else {
				A[i,j]=-(2.0/3.0)*dt*jacobian[i][j];
			}
		}
		b[i]=(4.0/3.0)*state[i] + (2.0/3.0)*dt*rate[i] - (1.0/3.0)*prev_state[i] - state_guess[i];
	}
	//solve system with bbla
	int[] pivot;
	pivot.length=state.length;
	LUDecomp!double(A,pivot);
	LUSolve!double(A,pivot,b,x);
	
	double[] update_guess1=state_guess.dup;
	update_guess1[] += x[];
	
	return update_guess1;
}

double[] Var_BDF2(double prev_dt, double dt, double[] prev_state, double[] state, double[] state_guess, double[] rate, double[][] jacobian) {
	//Solution takes the form Ax=b
	//[I-gamma*dt(dR/dS)](dS) = alpha*state + gamma*dt*rate - beta*prev_state - state_guess
	auto A = new Matrix!double(state.length,state.length);
	double[] x;
	x.length=state.length;
	double[] b;
	b.length=state.length;
	
	//define coefficients
	double w=(dt/prev_dt);
	double alpha=((1+w)^^2)/(1+2*w);
	double beta=(w^^2)/(1+2*w);
	double gamma=(1+w)/(1+2*w);
	
	//populate A and b
	for (int i;i<state.length;i++) {
		for (int j;j<state.length;j++) {
			if (i==j) {
				A[i,j]=1-gamma*dt*jacobian[i][j];
			} else {
				A[i,j]=-gamma*dt*jacobian[i][j];
			}
		}
		b[i]=alpha*state[i] + gamma*dt*rate[i] - beta*prev_state[i] - state_guess[i];
	}
	
	//solve system with bbla
	int[] pivot;
	pivot.length=state.length;
	LUDecomp!double(A,pivot);
	LUSolve!double(A,pivot,b,x);
	
	double[] update_guess=state_guess.dup;
	update_guess[] += x[];
	
	return update_guess;
}
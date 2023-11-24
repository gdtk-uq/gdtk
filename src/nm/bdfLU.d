/** 
* bdfLU.d
* 
* This module contains a set of implicit backwards differentiation methods.
* These methods are written to solve the matrix equations developed when 
* applying a newton-raphson step to a vector.
* These are an extension of the backward euler step, using a newton-raphson
* method to solve for the k+1 state.
*
* These methods call on Crout's LU decomposition method as specified in bbla.d
* to solve the set of matrix equations. These are structured with the intent of
* being called multiple times for any additional newton steps to improve the
* accuracy of each timestep.
*
* Contents:
* 	1. BDF1 - A single newton-raphson step using a first order backward
* 				differentiation formula.
* 	2. BDF2 - A single newton-raphson step using a second order backward
* 				differentiation formula.
* 	3. Var_BDF2 - A single newton-raphson step using a second order backward
* 				differentiation formula with a variable timestep.
*
* Author: Brad Semple
* Version: 2018-09-28, written
*/


module nm.bdfLU;

import std.stdio;
import std.math;
import std.conv;
import std.algorithm;
import std.exception;
import nm.bbla;
import ntypes.complex;
import nm.number;

//BDF1 or backwards euler formula
//applies newton-raphson method to solve
//requires current state vector, current rate vector, jacobian (dR/dS), and timestep
//returns new guess state vector for newton method

void BDF1(number[] update_vector, number[] RHSvector, Matrix!number Asolve, int[] pivot, double dt, number[] state, ref number[] state_guess, number[] rate, number[][] jacobian) {
	//Solution takes the form Ax=b
	//[I-dt(dR/dS)](dS) = state + dt*rate - state_guess
	//populate A and b
	for (int i;i<state.length;i++) {
		for (int j;j<state.length;j++) {
			if (i==j) {
				Asolve[i,j]=1-dt*jacobian[i][j];
			} else {
				Asolve[i,j]=-dt*jacobian[i][j];
			}
		}
		RHSvector[i]=state[i] + dt*rate[i] - state_guess[i];
	}
	
	LUDecomp!number(Asolve,pivot);
	LUSolve!number(Asolve,pivot,RHSvector,update_vector);

	foreach (int i; 0 .. to!int(update_vector.length)){
		state_guess[i] += update_vector[i];
	}
}

void BDF2(number[] update_vector, number[] RHSvector, Matrix!number Asolve, int[] pivot, double dt, number[] prev_state, number[] state, ref number[] state_guess, number[] rate, number[][] jacobian) {
	//Solution takes the form Ax=b
	//[I-(2/3)*dt(dR/dS)](dS) = (4/3)*state + (2/3)*dt*rate - (1/3)*prev_state - state_guess
	//populate A and b
	for (int i;i<state.length;i++) {
		for (int j;j<state.length;j++) {
			if (i==j) {
				Asolve[i,j]=1-(2.0/3.0)*dt*jacobian[i][j];
			} else {
				Asolve[i,j]=-(2.0/3.0)*dt*jacobian[i][j];
			}
		}
		RHSvector[i]=(4.0/3.0)*state[i] + (2.0/3.0)*dt*rate[i] - (1.0/3.0)*prev_state[i] - state_guess[i];
	}
	
	LUDecomp!number(Asolve,pivot);
	LUSolve!number(Asolve,pivot,RHSvector,update_vector);
	
	foreach (int i; 0 .. to!int(update_vector.length)){
		state_guess[i] += update_vector[i];
	}
}

void Var_BDF2(number[] update_vector, number[] RHSvector, Matrix!number Asolve, int[] pivot, double dt, double prev_dt, number[] prev_state, number[] state, ref number[] state_guess, number[] rate, number[][] jacobian) {
	//Solution takes the form Ax=b
	//[I-gamma*dt(dR/dS)](dS) = alpha*state + gamma*dt*rate - beta*prev_state - state_guess
	
	//define coefficients
	double w=(dt/prev_dt);
	double alpha=((1+w)^^2)/(1+2*w);
	double beta=(w^^2)/(1+2*w);
	double gamma=(1+w)/(1+2*w);
	
	//populate A and b
	for (int i;i<state.length;i++) {
		for (int j;j<state.length;j++) {
			if (i==j) {
				Asolve[i,j]=1-gamma*dt*jacobian[i][j];
			} else {
				Asolve[i,j]=-gamma*dt*jacobian[i][j];
			}
		}
		RHSvector[i]=alpha*state[i] + gamma*dt*rate[i] - beta*prev_state[i] - state_guess[i];
	}
	
	LUDecomp!number(Asolve,pivot);
	LUSolve!number(Asolve,pivot,RHSvector,update_vector);
	
	foreach (int i; 0 .. to!int(update_vector.length)){
		state_guess[i] += update_vector[i];
	}
}

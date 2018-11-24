/**
 * Author: Brad S.
 * Date: 2018-09-05
 * 
 * This numerical integrator provides a solution to coupled 
 * electronic species interactions through electron impact.
 * 
 * The current model considers binned states of Nitrogen and 
 * Oxygen from the ground state up to the 46th and 40th state
 * respectively.
 *
 */

module kinetics.electronic_state_solver;

import std.stdio;
import std.file;
import std.string;
import std.format;
import std.conv : to;
import std.math;
import std.algorithm.iteration : sum;

import nm.bbla;
import nm.bdfLU;
import nm.complex;
import nm.number;
import util.msg_service;


//define a list of global run conditions to take from file
number Te;
double endtime;
double dt;
double dj;
int newton_steps;
number Ne;

//define arrays for data storage for gas info and curve fitting parameters
double[8][47][47][2] full_rate_fit;
int total_species_number;

number[][] state;
number[][] guess_state;
number[][] step_state;
number[] state_vector;
number[] prev_state_vector;
number[] initial_vector;

//define linalg terms for BDF solver
number[] rate_vector;
number[][] jacobian;
Matrix!number Asolve;
number[] update_vector;
number[] RHSvector;
int[] pivot;

@nogc
void UpdateVectorFromState(ref number[] given_vector, number[][] given_state) 
{
    int j=0;
    foreach(number[] row; given_state) {
        foreach(int i, number elem; row) {
            given_vector[i+j] = elem;
        }
        j += row.length;
    }
}

@nogc
void UpdateStateFromVector(ref number[][] given_state, number[] given_vector)
{
    int Vector_Pos(int m, int n){
        int pos_sum=0;
        for(int i;i<m;i++) {
            pos_sum+=given_state[i].length;
        }
        pos_sum+=n;
        return pos_sum;
    }
    foreach(int m, number[] chem_spec;given_state) {
        foreach(int n, number elec_spec;chem_spec) {
            given_state[m][n]=given_vector[Vector_Pos(m,n)];
        }
    }
}

@nogc
void CopyStateToState(ref number[][] out_state, number[][] in_state, bool step=false, int stepn=0)
{
    int[2] Species_Position(int vector_position) {
        int chem_spec;
        while(vector_position>=out_state[chem_spec].length){
            vector_position-=out_state[chem_spec].length;
            chem_spec+=1;
        }
        return [chem_spec,vector_position];
    }

    foreach(int i, number[] row;in_state) {
        foreach (int j, number elem;row) {
            out_state[i][j] = elem;
        }
    }
    if (step) {
        int[2] step_pos=Species_Position(stepn);
        out_state[step_pos[0]][step_pos[1]] += dj;
    }
}

void Init(double[8][47][47][2] incoming_rate_fits, int[] species_numbers) 
{   
    full_rate_fit=incoming_rate_fits;
    //populate energy data and rate fitting coefficients

    debug {

        //set size of state arrays and linalg BDF solver arrays and vectors
        state.length=species_numbers.length;
        guess_state.length=species_numbers.length;
        step_state.length=species_numbers.length;

        foreach (int i, int indiv_data;species_numbers){
            total_species_number+=indiv_data+1;

            state[i].length=indiv_data+1;
            guess_state[i].length=indiv_data+1;
            step_state[i].length=indiv_data+1;
        }
        state_vector.length=total_species_number;
        prev_state_vector.length=total_species_number;
        initial_vector.length=total_species_number;
        rate_vector.length=total_species_number;
        Asolve = new Matrix!number(total_species_number,total_species_number);
        update_vector.length=total_species_number;
        RHSvector.length=total_species_number;
        pivot.length=total_species_number;

        //set size for linalg BDF solver arrays and vectors
        jacobian.length=total_species_number;
        foreach(int i,number[] row;jacobian){
            jacobian[i].length=total_species_number;
        }
    } else {
        throw new Error("PopulateRateFits is only available in debug flavour.");
    }
}

@nogc
number Rate(int ip, int gas, number[][] input_state=state, number electron_density=Ne) 
{
    //Returns the current rate of change of a grouped electronic species

    @nogc
    number Forward_Rate_Coef(int gas, int ip,int jp) {
        double A=full_rate_fit[gas][ip+1][jp+1][0];
        double n=full_rate_fit[gas][ip+1][jp+1][1];
        double E=full_rate_fit[gas][ip+1][jp+1][2];
        return A*((Te)^^n)*exp(-E/(Te));
    }

    @nogc
    number Equil_Const(int gas, int ip, int jp) {
        double G1=full_rate_fit[gas][ip+1][jp+1][3];
        double G2=full_rate_fit[gas][ip+1][jp+1][4];
        double G3=full_rate_fit[gas][ip+1][jp+1][5];
        double G4=full_rate_fit[gas][ip+1][jp+1][6];
        double G5=full_rate_fit[gas][ip+1][jp+1][7];
        number z=10000/(Te);
        return exp((G1/z) + G2 + G3*log(z) + G4*z + G5*(z^^2));
    }

    @nogc
    number Backward_Rate_Coef(int gas, int ip, int jp) {
        return Forward_Rate_Coef(gas,ip,jp)/Equil_Const(gas,ip,jp);
    }

    number group_rate=0;
	if (ip<(input_state[gas].length-1)) { //excited atomic species
		for (int jp;jp<(input_state[gas].length-1);jp++) {
			if (ip<jp) {
				group_rate+=Backward_Rate_Coef(gas,ip,jp)*input_state[gas][jp]*electron_density - Forward_Rate_Coef(gas,ip,jp)*input_state[gas][ip]*electron_density;
			} else if (ip>jp) {
				group_rate+=Forward_Rate_Coef(gas,jp,ip)*input_state[gas][jp]*electron_density - Backward_Rate_Coef(gas,jp,ip)*input_state[gas][ip]*electron_density;
			}
		}
		group_rate+=Backward_Rate_Coef(gas,ip,-1)*input_state[gas][$-1]*(electron_density^^2) - Forward_Rate_Coef(gas,ip,-1)*input_state[gas][ip]*electron_density;
	} else { //ionised species
		for (int jp;jp<(input_state[gas].length-1);jp++) {
			group_rate+=Forward_Rate_Coef(gas,jp,-1)*input_state[gas][jp]*electron_density - Backward_Rate_Coef(gas,jp,-1)*input_state[gas][ip]*(electron_density^^2);
		}
	}
	return group_rate;
}

@nogc
void Jacobian()
{
    @nogc
    bool Ion(int vector_position){
        int chem_spec;
        while(vector_position>=guess_state[chem_spec].length){
            vector_position-=guess_state[chem_spec].length;
            chem_spec+=1;
        }
        if (vector_position==guess_state[chem_spec].length-1) {
            return true;
        } else {
            return false;
        }
    }

    @nogc
    int[2] Species_Position(int vector_position) {
        int chem_spec;
        while(vector_position>=guess_state[chem_spec].length){
            vector_position-=guess_state[chem_spec].length;
            chem_spec+=1;
        }
        return [chem_spec,vector_position];
    }

    for(int n;n<total_species_number;n++) { //for every element of the state
        CopyStateToState(step_state,guess_state,true,n); //update step state
        number step_Ne=Ne;
        if (Ion(n)) {
            step_Ne+=dj;
        }
        for(int m;m<total_species_number;m++) { //for every element of the rate
            int[2] positions=Species_Position(m);
            int chem_spec=positions[0];
            int elec_spec=positions[1];
            jacobian[m][n]=(Rate(elec_spec,chem_spec,step_state,step_Ne)-Rate(elec_spec,chem_spec,guess_state,Ne))/dj;
        }
    }
}

@nogc
void Step() 
{   
    CopyStateToState(guess_state,state);
    UpdateVectorFromState(initial_vector,state);

    @nogc
    int VectorPosition(int gas, int ip)
    {
        int vectorposition;
        foreach(int i; 0 .. gas) {
            vectorposition+=state[i].length;
        }
        vectorposition+=ip;
        return vectorposition;
    }
    //step dt with BDF formula
    for (int n;n<newton_steps;n++) {
        foreach(int gas,number[] chem_spec;guess_state) {
            foreach(int ip,number elec_state;chem_spec){
                rate_vector[VectorPosition(gas,ip)]=Rate(ip, gas, guess_state);
            }
        }
        UpdateVectorFromState(state_vector,guess_state);
        number initial_ions=0;
        foreach(number[] chem_spec;guess_state) {
            initial_ions+=chem_spec[$-1];
        }
        Jacobian();
        debug {
            BDF2(update_vector, RHSvector, Asolve, pivot, dt, prev_state_vector, initial_vector, state_vector,rate_vector,jacobian);
            UpdateStateFromVector(guess_state,state_vector);
            number final_ions=0.0;
            foreach(number[] chem_spec;guess_state) {
                final_ions+=chem_spec[$-1];
            }
            Ne+=final_ions-initial_ions;
        } else {
            throw new Error("Var_BDF2 is only available in debug flavour.");
        }
    }
    CopyStateToState(state,guess_state);
    prev_state_vector[] = state_vector[];
}

@nogc
void First_Step()
{
    CopyStateToState(guess_state,state);
    UpdateVectorFromState(initial_vector,state);

    @nogc
    int VectorPosition(int gas, int ip) //A small function for give the position in the continuous vector, from a grouped species of a gas
    {
        int vectorposition;
        foreach(int i; 0 .. gas) {
            vectorposition+=state[i].length;
        }
        vectorposition+=ip;
        return vectorposition;
    }
    //step dt with BDF formula
    for (int n;n<newton_steps;n++) {
        foreach(int gas,number[] chem_spec;guess_state) {
            foreach(int ip,number elec_state;chem_spec){
                rate_vector[VectorPosition(gas,ip)]=Rate(ip, gas, guess_state);
            }
        }
        UpdateVectorFromState(state_vector,guess_state);
        number initial_ions=0;
        foreach(number[] chem_spec;guess_state) {
            initial_ions+=chem_spec[$-1];
        }
        Jacobian();
        debug {
            BDF1(update_vector, RHSvector, Asolve, pivot, dt/100.0,initial_vector,state_vector,rate_vector,jacobian);
            UpdateStateFromVector(guess_state,state_vector);
            number final_ions=0.0;
            foreach(number[] chem_spec;guess_state) {
                final_ions+=chem_spec[$-1];
            }
            Ne+=final_ions-initial_ions;
        } else {
            throw new Error("BDF1 is only available in debug flavour.");
        }
    }
    CopyStateToState(state,guess_state);
    prev_state_vector[] = state_vector[];

    for (int n;n<newton_steps;n++) {
        foreach(int gas,number[] chem_spec;guess_state) {
            foreach(int ip,number elec_state;chem_spec){
                rate_vector[VectorPosition(gas,ip)]=Rate(ip, gas, guess_state);
            }
        }
        UpdateVectorFromState(state_vector,guess_state);
        number initial_ions=0;
        foreach(number[] chem_spec;guess_state) {
            initial_ions+=chem_spec[$-1];
        }
        Jacobian();
        debug {
            Var_BDF2(update_vector, RHSvector, Asolve, pivot,(99.0/100.0)*dt, dt/100.0, prev_state_vector, initial_vector, state_vector,rate_vector,jacobian);
            UpdateStateFromVector(guess_state,state_vector);
            number final_ions=0.0;
            foreach(number[] chem_spec;guess_state) {
                final_ions+=chem_spec[$-1];
            }
            Ne+=final_ions-initial_ions;
        } else {
            throw new Error("Var_BDF2 is only available in debug flavour.");
        }
    }
    CopyStateToState(state,guess_state);
    prev_state_vector[] = state_vector[];
}

@nogc
int Electronic_Solve( number[] state_from_cfd, ref number[] state_to_cfd, number given_Te, double given_endtime, double dtChemSuggest)
{   
    debug{
        assert((given_Te>1500 && given_Te<30000.0), brokenPreCondition("given_Te"));
        assert((dtChemSuggest>1e-13 && dtChemSuggest<1e-5), brokenPreCondition("dtChemSuggest"));
        foreach(int i, number val;state_from_cfd) {
            assert ((!isNaN(val.re)),brokenPreCondition("state from cfd value"));
        }
    }

    // debug{
    //     writeln(state_from_cfd);
    //     writeln(dtChemSuggest);

    // }
    dt = dtChemSuggest;
    Ne=state_from_cfd[$-1];
    state_vector[]=state_from_cfd[0 .. $-1];
    UpdateStateFromVector(state,state_vector);
    Te=given_Te;
    endtime=given_endtime;
    dj=20;
    newton_steps=3;
    double t=0.0;
    if (sum(state_vector) < 1.0) {
        state_to_cfd[0 .. state_vector.length] = state_vector[];
        state_to_cfd[state_vector.length] = Ne;
        return 0;
    }
    First_Step();
    t+=dt;
    while (t<endtime) {
        Step();
        t+=dt;
        if ((endtime-t)<dt){
            dt = (endtime-t);
        }
        if (dt<1e-18){
            break;
        }
    }
    UpdateVectorFromState(state_vector,state);
    state_to_cfd[0 .. state_vector.length] = state_vector[];
    state_to_cfd[state_vector.length] = Ne;
    // debug{
    //     writeln(state_to_cfd);
    // }
    return 0;
}


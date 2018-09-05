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

import nm.bbla;
import nm.bdfLU;

//define a list of global run conditions to take from file
double Te;
double endtime;
double dt;
double dj;
int newton_steps;
double Ne;

//define arrays for data storage for gas info and curve fitting parameters
double[][][] full_grouped_data;
double[][string] N_rate_fit;
double[][string] O_rate_fit;
double[][string][] full_rate_fit;
double[] flow_code_state;
double[][] initial_state;
double[][] state;
double[] state_vector;

double[] Create_Vector(double[][] input_state) { //Creates vector of given state "input state"
    state_vector = [];
    foreach(double[] row;input_state) {
        state_vector~=row[];
    }
    return state_vector;
}

double[][] Create_State(double[] given_vector) { //Creates 2D state array from given vector
    int Vector_Pos(int m, int n){
        int pos_sum=0;
        for(int i;i<m;i++) {
            pos_sum+=state[i].length;
        }
        pos_sum+=n;
        return pos_sum;
    }
    double[][] output_state=Copy_State(state);
    foreach(int m, double[] chem_spec;state) {
        foreach(int n, double elec_spec;chem_spec) {
            output_state[m][n]=given_vector[Vector_Pos(m,n)];
        }
    }
    return output_state;
}

double[][] Copy_State(double[][] input_state) {//Creates copy of 2D state array from a given state
    double[][] output_data;
    output_data.length=input_state.length;
    foreach(int i, double[] row;input_state) {
        foreach (double elem;row) {
            output_data[i]~=elem;
        }
    }
    return output_data;
}

void Init() {
    //initialise electron density from last element of given state
    Ne=initial_state_vector[$-1];

    //create initial state and state
    initial_state.length=2;
    state.length=2;
    foreach(int i ; 0 .. 16) {
        if (i<8) {
            initial_state[0]~=initial_state_vector[i];
            state[0]~=initial_state_vector[i];
        } else {
            initial_state[1]~=initial_state_vector[i];
            state[1]~=initial_state_vector[i];
        }
    }
    //create state vector from state
    state_vector=Create_Vector(state);

    //initialise gas data from gas_data.d
    populate_rate_fits();
    full_grouped_data~=to!(double[][])(NI_grouped_data);
    full_grouped_data~=to!(double[][])(OI_grouped_data);
    full_rate_fit~=to!(double[][string])(N_rate_fit);
    full_rate_fit~=to!(double[][string])(O_rate_fit);
}

double Rate(int ip, int gas, double[][] input_state=state, double electron_density=Ne) {
    //Takes an energy level, gas and current state and returns the rate of change of grouped level ip\
    //if no state is specified, use current state
    //if no electron density is specified, use Ne
    double Forward_Rate_Coef(int gas, int ip,int jp) {
        double A=full_rate_fit[gas][format("%s,%s",ip+1,jp+1)][0];
        double n=full_rate_fit[gas][format("%s,%s",ip+1,jp+1)][1];
        double E=full_rate_fit[gas][format("%s,%s",ip+1,jp+1)][2];
        return A*((Te)^^n)*exp(-E/(Te));
    }

    double Equil_Const(int gas, int ip, int jp) {
        double G1=full_rate_fit[gas][format("%s,%s",ip+1,jp+1)][3];
        double G2=full_rate_fit[gas][format("%s,%s",ip+1,jp+1)][4];
        double G3=full_rate_fit[gas][format("%s,%s",ip+1,jp+1)][5];
        double G4=full_rate_fit[gas][format("%s,%s",ip+1,jp+1)][6];
        double G5=full_rate_fit[gas][format("%s,%s",ip+1,jp+1)][7];
        double z=10000/(Te);
        return exp((G1/z) + G2 + G3*log(z) + G4*z + G5*(z^^2));
    }

    double Backward_Rate_Coef(int gas, int ip, int jp) {
        return Forward_Rate_Coef(gas,ip,jp)/Equil_Const(gas,ip,jp);
    }

    double group_rate=0;
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

double[][] Jacobian(double[][] given_state){
    int size;
    foreach(double[] chem_spec;given_state) {
        size+=chem_spec.length;
    }
    double[][] jacobian;
    jacobian.length=size;
    foreach(int i,double[] row;jacobian){
        jacobian[i].length=size;
    }
    double[][] Create_Step_State(int vector_position) {
        double[][] output_state=Copy_State(given_state);
        int chem_spec;
        while (vector_position>=given_state[chem_spec].length) {
            chem_spec+=1;
            vector_position-=given_state[chem_spec].length;
        }
        int ip=vector_position;
        output_state[chem_spec][ip]+=dj;
        return output_state;
    }

    bool Ion(int vector_position){
        int chem_spec;
        while(vector_position>=given_state[chem_spec].length){
            chem_spec+=1;
            vector_position-=given_state[chem_spec].length;
        }
        if (vector_position==given_state[chem_spec].length-1) {
            return true;
        } else {
            return false;
        }
    }

    int[] Species_Position(int vector_position) {
        int chem_spec;
        while(vector_position>=given_state[chem_spec].length){
            chem_spec+=1;
            vector_position-=given_state[chem_spec].length;
        }
        return [chem_spec,vector_position];
    }

    for(int n;n<size;n++) { //for every element of the state
        double[][] step_state=Create_Step_State(n);
        double step_Ne=Ne;
        if (Ion(n)) {
            step_Ne+=dj;
        }
        for(int m;m<size;m++) { //for every element of the rate
            int[] positions=Species_Position(m);
            int chem_spec=positions[0];
            int elec_spec=positions[1];
            jacobian[m][n]=(Rate(elec_spec,chem_spec,step_state,step_Ne)-Rate(elec_spec,chem_spec,given_state,Ne))/dj;
        }
    }
    return jacobian;
}

void Step() {
    double[][] guess_state=Copy_State(state);
    double[][] initial_state=Copy_State(state);
    double[] initial_vector=Create_Vector(initial_state);

    //step dt with BDF formula
    for (int n;n<newton_steps;n++) {
        double[] rate_vector;
        foreach(int gas,double[] chem_spec;guess_state) {
            foreach(int ip,double elec_state;chem_spec){
                rate_vector~=Rate(ip, gas, guess_state);
            }
        }
        state_vector=Create_Vector(guess_state);
        double initial_ions=0;
        foreach(double[] chem_spec;guess_state) {
            initial_ions+=chem_spec[$-1];
        }
        writeln(rate_vector);
        writeln(state_vector);
        state_vector=BDF1(dt,initial_vector,state_vector,rate_vector,Jacobian(guess_state));
        guess_state=Create_State(state_vector);
        double final_ions=0.0;
        foreach(double[] chem_spec;guess_state) {
            final_ions+=chem_spec[$-1];
        }
        Ne+=final_ions-initial_ions;
    }
    state=Copy_State(guess_state);
}

void Electronic_Solve(in double[] state_from_cfd, out double[] state_to_cfd, double given_Te, double given_endtime)
{
    immutable double dt = 1.0e-10;
    initial_state_vector=state_from_cfd;
    Init();
    Te=given_Te;
    endtime=given_endtime;
    dj=20;
    newton_steps=3;

    double t=0.0;
    while (t<endtime) {
        Step();
        writeln(state_vector);
        t+=dt;
        if ((endtime-t)<dt){
            dt=(endtime-dt);
        }
    }

    // For consideration....
    // Avoid the "create vector"
    state_vector = Create_Vector(state);
    state_vector ~= Ne;
    state_to_cfd[] = state_vector[];
    return update_cfd_state;
}


// current grouping model for NI
// 1: 1
// 2: 2
// 3: 3
// 4: 4-6
// 5: 7-13
// 6: 14-21
// 7: 22-27
// 8: 28-46

immutable double[][] NI_grouped_data = [[1, 1, 1, 0, 4], 
                                        [2, 2, 2, 2.384, 10], 
                                        [3, 3, 3, 3.576, 6], 
                                        [4, 4, 6, 10.641, 30], 
                                        [5, 7, 13, 11.9508, 64], 
                                        [6, 14, 21, 12.9848, 110], 
                                        [7, 22, 27, 13.342, 64], 
                                        [8, 28, 46, 13.9876, 922]]; //grouped number, lower level, upper level, energy, degeneracy

// current grouping model for OI
// 1: 1
// 2: 2
// 3: 3
// 4: 4-6
// 5: 7-13
// 6: 14-21
// 7: 22-27
// 8: 28-40

immutable double[][] OI_grouped_data = [[1, 1, 1, 0, 9], 
                                        [2, 2, 2, 1.97, 5], 
                                        [3, 3, 3, 4.19, 1], 
                                        [4, 4, 6, 10.2345, 23], 
                                        [5, 7, 13, 12.0181, 81], 
                                        [6, 14, 21, 12.7522, 139], 
                                        [7, 22, 27, 13.0729, 128], 
                                        [8, 28, 40, 13.3392, 428]]; //grouped number, lower level, upper level, energy, degeneracy

void populate_rate_fits(){
    N_rate_fit = ["1,6":[3.62e-10, 0.354826, 152699, 0.000185762, 3.31603, 0.00146398, -15.0706, 0.00274539], 
                    "2,6":[2.26e-08, -0.00885559, 127549, 0.000185762, 2.39974, 0.00146398, -12.3041, 0.00274539], 
                    "3,6":[1.5e-08, 0.0505113, 113624, 0.000185762, 2.91057, 0.00146398, -10.9208, 0.00274539], 
                    "5,6":[4.39e-07, 0.0973882, 13334, 0.00543894, 0.601935, 0.0451428, -1.28467, -0.00796418], 
                    "4,6":[1.53e-07, 0.0607039, 29814.5, 0.0133102, 1.43796, 0.107322, -2.90595, -0.0105668], 
                    "2,7":[2.78e-09, 0.0932004, 126748, -0.00245775, 1.82911, -0.020368, -12.678, 0.00330493], 
                    "4,5":[1.09e-07, 0.121831, 15339.6, 0.00787124, 0.836029, 0.0621793, -1.62128, -0.00260263], 
                    "5,7":[6.09e-08, 0.115274, 17463.2, 0.00279543, 0.0313057, 0.0233108, -1.6586, -0.00740463], 
                    "2,5":[7.64e-08, -0.0975305, 113662, -0.00525318, 1.79781, -0.0436788, -11.0194, 0.0107096], 
                    "2,8":[4.76e-09, 0.19384, 132311, -0.00759934, 4.44976, -0.0594608, -13.3697, 0.015825], 
                    "3,8":[5.04e-09, 0.201329, 118460, -0.00759934, 4.96058, -0.0594608, -11.9864, 0.015825], 
                    "1,0":[5.69e-12, 0.721, 167900, 0.138661, 50.6526, -1.49173, -15.8033, -0.000702945], 
                    "2,3":[5.83e-06, -0.5364, 26677.9, -7.31e-14, -0.510826, -4.67e-13, -1.38326, -3.62e-14], 
                    "1,5":[2.74e-09, 0.204787, 136614, -0.00525318, 2.7141, -0.0436788, -13.786, 0.0107096], 
                    "4,7":[5.7e-08, 0.0630134, 33597.8, 0.0106667, 0.867334, 0.0854901,-3.27988, -0.0100073], 
                    "1,7":[3.44e-09, 0.0696977, 154330, -0.00245775, 2.7454, -0.020368, -15.4446, 0.00330493], 
                    "3,4":[1.61e-08, -0.0259361, 80756.4, -0.0131244, 1.4726, -0.105858, -8.0149, 0.0133122], 
                    "6,8":[4.06e-07, 0.185587, 9283.29, -0.0077851, 2.05002, -0.0609248, -1.06557, 0.0130796], 
                    "4,0":[1.15e-09, 0.48233, 48055.7, 0.151785, 48.7745, -1.38587, -3.63864, -0.0140151], 
                    "8,0":[2.32e-05, -0.101767, 9610.63, 0.14626, 45.2866, -1.43227, 0.332881, -0.0165279], 
                    "7,0":[2.77e-07, 0.173263, 15499.1, 0.141119, 47.9072, -1.47136, -0.358756, -0.00400788], 
                    "1,3":[3.11e-08, -0.204736, 44700.5, 4.56e-14, 0.405465, 3.11e-13, -4.14978, 2.87e-14], 
                    "5,0":[1.29e-08, 0.33318, 32293.8, 0.143914, 47.9385, -1.44805, -2.01735, -0.0114125], 
                    "1,2":[5.77e-08, -0.160269, 40116.9, 4.81e-14, 0.916291, 3.45e-13, -2.76652, 3.63e-14], 
                    "3,0":[1.62e-11,0.6691, 126600, 0.138661, 50.2471, -1.49173, -11.6535, -0.000702945], 
                    "5,8":[1.18e-07, 0.158589, 22725.6, -0.00234616, 2.65195, -0.015782, -2.35023, 0.00511541], 
                    "1,4":[7.99e-07, -0.320227, 117020, -0.0131244, 1.87807, -0.105858, -12.1647, 0.0133122], 
                    "4,8":[9.58e-08, 0.110706, 39225.5, 0.00552507, 3.48798, 0.0463973, -3.97152, 0.00251278], 
                    "2,4":[1.7e-08, -0.0430612, 94407.3, -0.0131244, 0.961777, -0.105858, -9.39816, 0.0133122], 
                    "7,8":[2.98e-06, 0.146262, 4920.13, -0.00514159, 2.62065, -0.0390928, -0.691636, 0.01252], 
                    "1,8":[4.54e-09, 0.178276, 160012, -0.00759934, 5.36605, -0.0594608, -16.1362, 0.015825], 
                    "3,7":[2.47e-09, 0.108029, 112984, -0.00245775, 2.33994, -0.020368, -11.2948, 0.00330493], 
                    "6,0":[5.99e-08, 0.269394, 19115, 0.138475, 47.3366, -1.49319, -0.732687, -0.00344833], 
                    "6,7":[1.35e-06, 0.144078, 3705.03, -0.00264351, -0.570629, -0.0218319, -0.373932, 0.000559542], 
                    "2,0":[7.21e-11, 0.5487, 140100, 0.138661, 49.7363, -1.49173, -13.0368, -0.000702945], 
                    "3,5":[3.77e-08, -0.030555, 97607.4, -0.00525318, 2.30863, -0.0436788, -9.63618, 0.0107096]];

    O_rate_fit = ["1,6":[3.76e-09, 0.0318765, 145547, 0.000989389, 2.74767, 0.00801344, -14.8125, 0.00732514], 
                    "2,6":[7.79e-09, 0.0519045, 123441, 0.000989389, 3.33546, 0.00801344, -12.5265, 0.00732514], 
                    "3,6":[3.32e-09, 0.145889, 98835.9, 0.000989389, 4.9449, 0.00801344, -9.95025, 0.00732514], 
                    "5,6":[4.07e-05, -0.210725, 17719.5, 0.0755789, 1.26922,0.589621, -1.75675, -0.0261783], 
                    "4,6":[2.24e-06, -0.165891, 40660.9, 0.0679141, 2.69173, 0.614996, -4.25644, 0.00428068], 
                    "2,7":[1.62e-08, -0.0285564, 128252, 0.000856145, 3.25118, 0.00677486, -12.8958, 0.00416073], 
                    "4,5":[7.73e-08, 0.15301, 23562.4, -0.00766474, 1.42251, 0.0253754, -2.49969, 0.030459], 
                    "5,7":[1.18e-05, -0.216561, 22015.2, 0.0754456, 1.18494, 0.588382, -2.12606, -0.0293428], 
                    "2,5":[9.01e-11, 0.409249, 110013, -0.0745895, 2.06624, -0.581607, -10.7697, 0.0335035], 
                    "2,8":[5.99e-09, 0.169608, 131515, -0.000582163, 4.44387, -0.00459808, -13.1857, 0.00471003], 
                    "3,8":[8.57e-09, 0.168803, 105836, -0.000582163, 6.05331, -0.00459808, -10.6095, 0.00471003], 
                    "1,0":[5.81e-11, 0.393085, 160067, 0.107316, 48.7824, -1.98736, -15.4748, -0.0199158], 
                    "2,3":[1.47e-09, 0.0156301, 26862.8, -2.33e-13, -1.60944, -1.48e-12, -2.5762, -1.18e-13], 
                    "1,5":[1.21e-09, 0.125162, 126957, -0.0745895, 1.47845, -0.581607, -13.0558, 0.0335035], 
                    "4,7":[1.23e-06, -0.17678, 44883.4, 0.0677809, 2.60745, 0.613758, -4.62575, 0.00111627], 
                    "1,7":[6.86e-10, 0.259424, 150830, 0.000856145, 2.66339, 0.00677486, -15.1819, 0.00416073], 
                    "3,4":[2.52e-09, 0.203356, 58130.8, -0.0669248, 2.25316, -0.606983, -5.69381, 0.00304446], 
                    "6,8":[9.11e-07, 0.150587, 6296.1, -0.00157155, 1.10842, -0.0126115, -0.65928, -0.00261511], 
                    "4,0":[9.69e-10, 0.230059, 53489.9, 0.17424, 48.7265, -1.38038, -4.91868, -0.0229603], 
                    "8,0":[2.14e-07, -0.140952, 5794.8, 0.107898, 44.9263, -1.98276, -0.00295947, -0.0246258], 
                    "7,0":[2.99e-08, 0.0117041, 9180.3, 0.10646, 46.119, -1.99413, -0.292929, -0.0240765], 
                    "1,3":[6.29e-12, 0.412098, 49014.7, -1.36e-13, -2.19722, -9.5e-13, -4.86229, -7.8e-14], 
                    "5,0":[4.08e-08, -0.0555187, 31478.7, 0.181905, 47.304, -1.40575, -2.41899, -0.0534193], 
                    "1,2":[1.57e-10, 0.312407, 23306.8, -4.83e-13, -0.587787, -3.2e-12, -2.28609, -2.95e-13], 
                    "3,0":[3.3e-13, 0.829595, 110776, 0.107316, 50.9796, -1.98736, -10.6125, -0.0199158], 
                    "5,8":[1.37e-05, -0.205288, 25151, 0.0740073, 2.37764, 0.577009, -2.41603, -0.0287934], 
                    "1,4":[1.95e-08, -0.146963, 106345, -0.0669248, 0.0559377, -0.606983, -10.5561, 0.00304446], 
                    "4,8":[2.27e-06, -0.160608, 48023.7, 0.0663426, 3.80015, 0.602385, -4.91572, 0.00166557], 
                    "2,4":[1.09e-09, 0.174335, 83260.7, -0.0669248, 0.643724, -0.606983, -8.27002, 0.00304446], 
                    "7,8":[6.69e-06, 0.15114, 1624.13, -0.00143831, 1.1927, -0.0113729, -0.289969, 0.000549307], 
                    "1,8":[4.62e-09, 0.17016, 154292, -0.000582163, 3.85609, -0.00459808, -15.4718, 0.00471003], 
                    "3,7":[4.27e-08, -0.0695154, 102549, 0.000856145, 4.86061, 0.00677486, -10.3196, 0.00416073], 
                    "6,0":[6.85e-09, 0.124579, 13163.1, 0.106326, 46.0347, -1.99537, -0.662239, -0.0272409], 
                    "6,7":[2.07e-06, 0.14048, 2752.66, -0.000133243, -0.0842826, -0.00123858, -0.36931, -0.00316441], 
                    "2,0":[1.67e-13, 0.867141, 136498, 0.107316, 49.3702, -1.98736, -13.1887, -0.0199158], 
                    "3,5":[1.67e-10,0.420516, 84171.4, -0.0745895, 3.67567, -0.581607, -8.1935, 0.0335035]];
}

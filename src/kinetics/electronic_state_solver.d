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
import nm.complex;
import nm.number;


//define a list of global run conditions to take from file
number Te;
double endtime;
double dt;
double dj;
int newton_steps;
number Ne;

//define arrays for data storage for gas info and curve fitting parameters
double[][][] full_grouped_data;
double[8][9][9][2] full_rate_fit;
int total_species_number;

number[][] state;
number[][] guess_state;
number[][] step_state;
number[] state_vector;
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
            chem_spec+=1;
            vector_position-=out_state[chem_spec].length;
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

void Init() 
{   
    //populate energy data and rate fitting coefficients
    full_grouped_data~=to!(double[][])(NI_grouped_data);
    full_grouped_data~=to!(double[][])(OI_grouped_data);
    PopulateRateFits();

    //set size of state arrays and linalg BDF solver arrays and vectors
    state.length=full_grouped_data.length;
    guess_state.length=full_grouped_data.length;
    step_state.length=full_grouped_data.length;

    foreach (int i, double[][] indiv_data;full_grouped_data){
        total_species_number+=indiv_data.length+1;

        state[i].length=indiv_data.length+1;
        guess_state[i].length=indiv_data.length+1;
        step_state[i].length=indiv_data.length+1;
    }
    state_vector.length=total_species_number;
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

}

@nogc
number Rate(int ip, int gas, number[][] input_state=state, number electron_density=Ne) 
{
    //Takes an energy level, gas and current state and returns the rate of change of grouped level ip\
    //if no state is specified, use current state
    //if no electron density is specified, use Ne

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
            chem_spec+=1;
            vector_position-=guess_state[chem_spec].length;
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
            chem_spec+=1;
            vector_position-=guess_state[chem_spec].length;
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
            BDF1(update_vector, RHSvector, Asolve, pivot, dt,initial_vector,state_vector,rate_vector,jacobian);
        } else {
            throw new Error("BDF1 is only available in debug flavour.");
        }
        UpdateStateFromVector(guess_state,state_vector);
        number final_ions=0.0;
        foreach(number[] chem_spec;guess_state) {
            final_ions+=chem_spec[$-1];
        }
        Ne+=final_ions-initial_ions;
    }
    CopyStateToState(state,guess_state);
}

@nogc
void Electronic_Solve( number[] state_from_cfd, ref number[] state_to_cfd, number given_Te, double given_endtime)
{
    dt = 1.0e-10;
    Ne=state_from_cfd[$-1];
    debug { writeln("state_vector: ",state_vector); }
    state_vector[]=state_from_cfd[0 .. $-1];
    UpdateStateFromVector(state,state_from_cfd);
    UpdateVectorFromState(state_vector,state);
    Te=given_Te;
    endtime=given_endtime;
    dj=20;
    newton_steps=3;
    double t=0.0;
    while (t<endtime) {
        Step();
        t+=dt;
        if ((endtime-t)<dt){
            dt=(endtime-dt);
        }
    }
    UpdateVectorFromState(state_vector,state);
    state_to_cfd[0 .. 18] = state_vector[];
    state_to_cfd[18] = Ne;
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

double[][] N_rate_fit = [[1,2,5.77E-08,-0.160268851,40116.89956,4.81E-14,0.916290732,3.45E-13,-2.766518071,3.63E-14],
                        [1,3,3.11E-08,-0.204736499,44700.53167,4.56E-14,0.405465108,3.11E-13,-4.149777107,2.87E-14],
                        [1,4,7.99E-07,-0.320227397,117019.9673,-0.013124412,1.878067796,-0.105858089,-12.16467391,0.013312199],
                        [1,5,2.74E-09,0.204786881,136614.3108,-0.005253176,2.714096489,-0.043678787,-13.78595741,0.010709565],
                        [1,6,3.62E-10,0.354825932,152698.7668,0.000185762,3.316031078,0.001463981,-15.07062382,0.002745391],
                        [1,7,3.44E-09,0.069697658,154329.7031,-0.002457748,2.745402204,-0.020367958,-15.44455544,0.003304932],
                        [1,8,4.54E-09,0.178275569,160011.6202,-0.007599338,5.366048321,-0.059460805,-16.13619155,0.015824978],
                        [2,3,5.83E-06,-0.536399907,26677.92703,-7.31E-14,-0.510825624,-4.67E-13,-1.383259036,-3.62E-14],
                        [2,4,1.70E-08,-0.043061242,94407.31075,-0.013124412,0.961777064,-0.105858089,-9.398155841,0.013312199],
                        [2,5,7.64E-08,-0.097530509,113662.3433,-0.005253176,1.797805757,-0.043678787,-11.01943934,0.010709565],
                        [2,6,2.26E-08,-0.008855589,127548.6113,0.000185762,2.399740346,0.001463981,-12.30410575,0.002745391],
                        [2,7,2.78E-09,0.093200448,126747.8289,-0.002457748,1.829111472,-0.020367958,-12.67803737,0.003304932],
                        [2,8,4.76E-09,0.193839649,132311.405,-0.007599338,4.449757589,-0.059460805,-13.36967348,0.015824978],
                        [3,4,1.61E-08,-0.02593611,80756.38125,-0.013124412,1.472602688,-0.105858089,-8.014896806,0.013312199],
                        [3,5,3.77E-08,-0.030554973,97607.36288,-0.005253176,2.308631381,-0.043678787,-9.636180306,0.010709565],
                        [3,6,1.50E-08,0.050511256,113623.7405,0.000185762,2.91056597,0.001463981,-10.92084671,0.002745391],
                        [3,7,2.47E-09,0.108028805,112984.1379,-0.002457748,2.339937096,-0.020367958,-11.29477833,0.003304932],
                        [3,8,5.04E-09,0.201329167,118460.4623,-0.007599338,4.960583213,-0.059460805,-11.98641444,0.015824978],
                        [4,5,1.09E-07,0.121830963,15339.60001,0.007871236,0.836028693,0.062179303,-1.6212835,-0.002602633],
                        [4,6,1.53E-07,0.060703853,29814.49084,0.013310174,1.437963282,0.107322071,-2.905949907,-0.010566808],
                        [4,7,5.70E-08,0.063013412,33597.79931,0.010666664,0.867334408,0.085490132,-3.279881529,-0.010007266],
                        [4,8,9.58E-08,0.110705527,39225.53684,0.005525074,3.487980525,0.046397285,-3.971517636,0.002512779],
                        [5,6,4.39E-07,0.097388243,13333.98025,0.005438938,0.601934589,0.045142768,-1.284666408,-0.007964175],
                        [5,7,6.09E-08,0.115274455,17463.16279,0.002795427,0.031305715,0.023310829,-1.658598029,-0.007404633],
                        [5,8,1.18E-07,0.158588783,22725.64975,-0.002346162,2.651951832,-0.015782018,-2.350234136,0.005115412],
                        [6,7,1.35E-06,0.144078277,3705.026149,-0.00264351,-0.570628874,-0.021831939,-0.373931622,0.000559542],
                        [6,8,4.06E-07,0.185587132,9283.292505,-0.0077851,2.050017243,-0.060924786,-1.065567729,0.013079587],
                        [7,8,2.98E-06,0.146262041,4920.133995,-0.00514159,2.620646117,-0.039092847,-0.691636107,0.012520046],
                        [1,0,5.69E-12,0.721,167900,0.138660786,50.6526134,-1.491728735,-15.80331096,-0.000702945],
                        [2,0,7.21E-11,0.5487,140100,0.138660786,49.73632267,-1.491728735,-13.03679289,-0.000702945],
                        [3,0,1.62E-11,0.6691,126600,0.138660786,50.24714829,-1.491728735,-11.65353386,-0.000702945],
                        [4,0,1.15E-09,0.482330444,48055.71996,0.151785197,48.77454561,-1.385870645,-3.638637053,-0.014015144],
                        [5,0,1.29E-08,0.333179859,32293.77081,0.143913961,47.93851691,-1.448049948,-2.017353553,-0.01141251],
                        [6,0,5.99E-08,0.269393946,19114.96237,0.138475024,47.33658232,-1.493192716,-0.732687145,-0.003448335],
                        [7,0,2.77E-07,0.173263358,15499.05671,0.141118534,47.9072112,-1.471360777,-0.358755523,-0.004007877],
                        [8,0,2.32E-05,-0.101766537,9610.631274,0.146260124,45.28656508,-1.43226793,0.332880584,-0.016527923]];

double[][] O_rate_fit = [[1,2,1.57E-10,0.312407211,23306.79906,-4.83E-13,-0.587786665,-3.20E-12,-2.286090856,-2.95E-13],
                        [1,3,6.29E-12,0.412098209,49014.74116,-1.36E-13,-2.197224577,-9.50E-13,-4.862294764,-7.80E-14],
                        [1,4,1.95E-08,-0.146962917,106345.2252,-0.066924752,0.055937729,-0.606982736,-10.55610657,0.003044461],
                        [1,5,1.21E-09,0.125161601,126956.639,-0.074589496,1.478449731,-0.581607386,-13.0557963,0.03350348],
                        [1,6,3.76E-09,0.031876532,145547.1461,0.000989389,2.747671171,0.00801344,-14.81254879,0.007325143],
                        [1,7,6.86E-10,0.259423861,150829.9274,0.000856145,2.663388556,0.006774861,-15.18185892,0.004160728],
                        [1,8,4.62E-09,0.170160191,154291.8281,-0.000582163,3.856087941,-0.004598076,-15.47182833,0.004710035],
                        [2,3,1.47E-09,0.015630071,26862.82121,-2.33E-13,-1.609437912,-1.48E-12,-2.576203909,-1.18E-13],
                        [2,4,1.09E-09,0.174335004,83260.69777,-0.066924752,0.643724394,-0.606982736,-8.270015717,0.003044461],
                        [2,5,9.01E-11,0.409248975,110012.7133,-0.074589496,2.066236396,-0.581607386,-10.76970544,0.03350348],
                        [2,6,7.79E-09,0.051904472,123440.745,0.000989389,3.335457836,0.00801344,-12.52645793,0.007325143],
                        [2,7,1.62E-08,-0.028556434,128252.4674,0.000856145,3.251175221,0.006774861,-12.89576806,0.004160728],
                        [2,8,5.99E-09,0.16960817,131515.1798,-0.000582163,4.443874606,-0.004598076,-13.18573747,0.004710035],
                        [3,4,2.52E-09,0.203355701,58130.78602,-0.066924752,2.253162307,-0.606982736,-5.693811808,0.003044461],
                        [3,5,1.67E-10,0.420515877,84171.44624,-0.074589496,3.675674308,-0.581607386,-8.193501531,0.03350348],
                        [3,6,3.32E-09,0.145888866,98835.94802,0.000989389,4.944895748,0.00801344,-9.950254023,0.007325143],
                        [3,7,4.27E-08,-0.069515374,102549.1392,0.000856145,4.860613134,0.006774861,-10.31956415,0.004160728],
                        [3,8,8.57E-09,0.168802758,105835.8906,-0.000582163,6.053312518,-0.004598076,-10.60953356,0.004710035],
                        [4,5,7.73E-08,0.15300974,23562.37452,-0.007664745,1.422512002,0.02537535,-2.499689723,0.030459019],
                        [4,6,2.24E-06,-0.16589096,40660.91695,0.06791414,2.691733442,0.614996176,-4.256442215,0.004280682],
                        [4,7,1.23E-06,-0.176779503,44883.4481,0.067780897,2.607450827,0.613757597,-4.625752344,0.001116267],
                        [4,8,2.27E-06,-0.16060752,48023.67295,0.066342588,3.800150212,0.602384661,-4.915721753,0.001665574],
                        [5,6,4.07E-05,-0.210724848,17719.52936,0.075578885,1.26922144,0.589620827,-1.756752492,-0.026178337],
                        [5,7,1.18E-05,-0.216561336,22015.2422,0.075445642,1.184938826,0.588382247,-2.126062621,-0.029342752],
                        [5,8,1.37E-05,-0.205288022,25151.03992,0.074007333,2.37763821,0.577009311,-2.416032031,-0.028793446],
                        [6,7,2.07E-06,0.140480078,2752.659364,-0.000133243,-0.084282614,-0.00123858,-0.369310128,-0.003164415],
                        [6,8,9.11E-07,0.150587181,6296.097376,-0.001571552,1.10841677,-0.012611516,-0.659279538,-0.002615108],
                        [7,8,6.69E-06,0.151139847,1624.127968,-0.001438308,1.192699384,-0.011372936,-0.28996941,0.000549307],
                        [1,0,5.81E-11,0.393085133,160066.8565,0.10731567,48.7824081,-1.987359243,-15.4747878,-0.019915794],
                        [2,0,1.67E-13,0.867140835,136498.0925,0.10731567,49.37019476,-1.987359243,-13.18869694,-0.019915794],
                        [3,0,3.30E-13,0.829595302,110775.8737,0.10731567,50.97963268,-1.987359243,-10.61249303,-0.019915794],
                        [4,0,9.69E-10,0.230058895,53489.85742,0.174240422,48.72647037,-1.380376507,-4.918681226,-0.022960255],
                        [5,0,4.08E-08,-0.055518737,31478.65713,0.181905167,47.30395837,-1.405751857,-2.418991504,-0.053419275],
                        [6,0,6.85E-09,0.124578918,13163.09466,0.106326281,46.03473693,-1.995372684,-0.662239011,-0.027240937],
                        [7,0,2.99E-08,0.011704109,9180.303463,0.106459525,46.11901954,-1.994134104,-0.292928883,-0.024076522],
                        [8,0,2.14E-07,-0.140952074,5794.796005,0.107897833,44.92632016,-1.982761168,-0.002959473,-0.024625829]];
void PopulateRateFits()
{
    foreach (double[] row;N_rate_fit){
        full_rate_fit[0][to!int(row[0])][to!int(row[1])] = row[2 .. $];
    }
    foreach (double[] row;O_rate_fit){
        full_rate_fit[1][to!int(row[0])][to!int(row[1])] = row[2 .. $];
    }
}

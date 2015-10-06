/**
 * file: unitproc.d
 * location: ./src/moc
 * author: Momar Hughes
 * description: Basic unit processes.
 * Each function computes the data for a new node based 
 * on information from other nodes.
 * The functions are computationally intensive and access 
 * the internals of the node data structure directly.
 * version: 	
 *		20 Mar 2015: initial port
 *		12 May 2015: use existing classes in cfcfd3
 *		22 Aug 2015: non-isentropic unit processes added
 *		16 Sep 2015: updated to accept new NodeData 
 *					 and FlowState classes
 */
module unitproc;

import std.stdio, std.math, std.conv, std.string;
import kernel;
import gasdynamic;
import wall;
import geom;
import gpath;

enum int MAX_ITERATION = 15; //maximum iteration count
enum int CORRECTOR_ITERATION = 4; //maximum iteration of corrector
enum double POS_TOL = 1.0e-5; //tolerance on positions
enum double PRES_TOL = 1.0; //tolerance on positions
//
int ISENTROPIC_FLAG = YES; // isentropic processes

void SetIsentropicFlag(int value)
{
    if(value == YES || value == NO){
	ISENTROPIC_FLAG = value;
    } else {
	throw new Error(text("Invalid flag, set INSENTROPIC FLAG as YES or NO"));
    } // end if
} // end SetIsentropicFlag()

int GetIsentropicFlag()
{
    return ISENTROPIC_FLAG;
} // end GetIsentropicFlag()

//************************************************************//
// Homentropic unit processes                                 //
// all from Zucrow M, Hoffman J, 1985, Gas Dynamics:          //
// Multidimensional Flow, Volume II, John Wiley and Sons      //
//************************************************************//

/**
 * Calculate an interior point from two initial points
 * Input:
 *   node1: index of initial point along C- characteristic
 *   node2: index of initial point along C+ characteristic
 *   node4: index of solution point (may have a value of -1)
 * If -1 is specified as the index for node4, a new node will be
 * created for the solution point.
 *   fs: flow state to work in
 * Output: 
 * Returns the index of the solution point 
 */
int InteriorNode(int node1,int node2,int node4,int fs=-1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    //
    if(GetIsentropicFlag==YES){
	return InteriorNode_0(node1,node2,node4,fs);
    } // end if
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //first, assign the Node data
    if(node1 == node2){
	throw new Error(text("InteriorNode: node1 & node2 have same id: ",node1));
    } // end if
    NodeData n1 = GetNodeData(node1);
    if (n1 is null){
	throw new Error(text("InteriorNode: node1 doesn't exist, id: ",node1));
    } // end if
    NodeData n2 = GetNodeData(node2);
    if (n2 is null){
	throw new Error(text("InteriorNode: node2 doesn't exist, id: ",node2));
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("StepStreamNode: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if 
    int node3 = InsertNode(node1,node2,-1,0.5);
    NodeData n3 = GetNodeData(node3);
    if (n3 is null){
	throw new Error(text("InteriorNode: could not create node3, id: ",node3));
    } // end if
    if (node4 == -1){
	node4 = CreateNode(-1);
    } // end if
    // define gas properties
    double mu1 = MachAngle(n1.Mach);
    double mu2 = MachAngle(n2.Mach);
    double mu3 = MachAngle(n3.Mach);
    // predictor variables
    double lambdaPlus = tan(n2.theta+mu2);
    double QPlus = sqrt(n2.Mach^^2-1.)/n2.rho/(n2.V^^2);
    double SPlus = axiflag*sin(n2.theta)/n2.y/n2.Mach/cos(n2.theta+mu2); //TODO fix for when y=0
    double lambdaMinus = tan(n1.theta-mu1);
    double QMinus = sqrt(n1.Mach^^2-1.)/n1.rho/(n1.V^^2);
    double SMinus = axiflag*sin(n1.theta)/n1.y/n1.Mach/cos(n1.theta-mu1);
    //location of point 4, for predictor
    n4.x = (n1.y-n2.y-lambdaMinus*n1.x+lambdaPlus*n2.x)/(lambdaPlus-lambdaMinus);
    n4.y = n2.y-lambdaPlus*(n2.x-n4.x);
    //
    double TPlus = -SPlus*(n4.x-n2.x)+QPlus*n2.P+n2.theta;
    double TMinus = -SMinus*(n4.x-n1.x)+QMinus*n1.P-n1.theta;
    //location of point 3, for predictor
    double lambda0 = tan(0.5*(n1.theta+n2.theta));
    double lambda12 = (n1.y-n2.y)/(n1.x-n2.x);
    //
    Vector3 xy3_old;
    double alpha;
    int iteration_count=0;
    double delta_position=0.0;
    do {
	++iteration_count;
	xy3_old = n3.pos;
	n3.x = (n2.y-n4.y-lambda12*n2.x+lambda0*n4.x)/(lambda0-lambda12);
	n3.y = n4.y-lambda0*(n4.x-n3.x);
	alpha = (n1.y-n2.y)/(n3.y-n2.y);
	n3.theta = n2.theta+(n1.theta-n2.theta)/alpha;
	lambda0 = tan(n3.theta);
	delta_position = abs(xy3_old-n3.pos);
    } while (iteration_count<MAX_ITERATION && delta_position>POS_TOL);
    //remaining flow properties of point 3, for predictor
    n3.P = n2.P+(n1.P-n2.P)/alpha;
    n3.rho = n2.rho+(n1.rho-n2.rho)/alpha;
    n3.T = n2.T+(n1.T-n2.T)/alpha;
    n3.V = n2.V+(n1.V-n2.V)/alpha;
    //flow properties of point 4, for predictor
    double R0 = n3.rho*n3.V;
    double A0 = n3.a^^2;
    double Tee01 = R0*n3.V + n3.P;
    double Tee02 = n3.P - A0*n3.rho;
    //
    n4.P = (TPlus+TMinus)/(QPlus+QMinus);
    n4.theta = TPlus-QPlus*n4.P;
    n4.V = (Tee01-n4.P)/R0;
    n4.rho = (n4.P-Tee02)/A0;
    n4.T = n4.P/R/n4.rho;
    //solve for corrector
    int count=0;
    double P4_old,delta_p_corrector=0.0;
    do{
	count++;
	P4_old = n4.P;
	//coefficients, for corrector
	double pPlus = 0.5*(n2.P+n4.P);
	double thetaPlus = 0.5*(n2.theta+n4.theta);
	double VPlus = 0.5*(n2.V+n4.V);
	double rhoPlus = 0.5*(n2.rho+n4.rho);
	double yPlus = 0.5*(n2.y+n4.y);
	double aPlus = sqrt(g*pPlus/rhoPlus);
	double MPlus = VPlus/aPlus;
	double muPlus = MachAngle(MPlus);
	lambdaPlus = tan(thetaPlus+muPlus);
	QPlus = sqrt(MPlus^^2-1.)/rhoPlus/VPlus^^2;
	SPlus = axiflag*sin(thetaPlus)/yPlus/MPlus/cos(thetaPlus+muPlus);
	//
	double pMinus = 0.5*(n1.P+n4.P);
	double thetaMinus = 0.5*(n1.theta+n4.theta);
	double VMinus = 0.5*(n1.V+n4.V);
	double rhoMinus = 0.5*(n1.rho+n4.rho);
	double yMinus = 0.5*(n1.y+n4.y);
	double aMinus = sqrt(g*pMinus/rhoMinus);
	double MMinus = VMinus/aMinus;
	double muMinus = MachAngle(MMinus);
	lambdaMinus = tan(thetaMinus-muMinus);
	QMinus = sqrt(MMinus^^2-1.)/rhoMinus/VMinus^^2;
	SMinus = axiflag*sin(thetaMinus)/yMinus/MMinus/cos(thetaMinus-muMinus);
	//location of point 4, for corrector
	n4.x = (n1.y-n2.y-lambdaMinus*n1.x+lambdaPlus*n2.x)/(lambdaPlus-lambdaMinus);
	n4.y = n2.y-lambdaPlus*(n2.x-n4.x);
	//
	TPlus = -SPlus*(n4.x-n2.x)+QPlus*n2.P+n2.theta;
	TMinus = -SMinus*(n4.x-n1.x)+QMinus*n1.P-n1.theta;
	//location of point 3, for corrector
	lambda0 = tan(0.5*(n3.theta+n4.theta));
	//
	iteration_count=0;
	delta_position=0.0;
	do {
	    ++iteration_count;
	    Vector3 xy3 = n3.pos;
	    n3.x = (n2.y-n4.y-lambda12*n2.x+lambda0*n4.x)/(lambda0-lambda12);
	    n3.y = n4.y-lambda0*(n4.x-n3.x);
	    alpha = (n1.y-n2.y)/(n3.y-n2.y);
	    n3.theta = n2.theta+(n1.theta-n2.theta)/alpha;
	    lambda0 = tan(n3.theta);
	    delta_position = abs(xy3-n3.pos);
	} while (iteration_count<MAX_ITERATION && delta_position>POS_TOL);
	//remaining flow properties of point 3, for corrector
	n3.P = n2.P+(n1.P-n2.P)/alpha;
	n3.rho = n2.rho+(n1.rho-n2.rho)/alpha;
	n3.T = n2.T+(n1.T-n2.T)/alpha;
	n3.V = n2.V+(n1.V-n2.V)/alpha;
	//flow properties of point 4, for corrector
	double p3 = 0.5*(n3.P+n4.P);
	double rho3 = 0.5*(n3.rho+n4.rho);
	double a3 = sqrt(g*p3/rho3);
	double V3 = 0.5*(n3.V+n4.V);
	R0 = rho3*V3;
	A0 = a3^^2;
	Tee01 = R0*V3 + p3;
	Tee02 = p3 - A0*rho3;
	//
	n4.P = (TPlus+TMinus)/(QPlus+QMinus);
	n4.theta = TPlus-QPlus*n4.P;
	n4.V = (Tee01-n4.P)/R0;
	n4.rho = (n4.P-Tee02)/A0;
	n4.T = n4.P/n4.rho/R;
	//
	delta_p_corrector = fabs(n4.P - P4_old);
    }while (count<CORRECTOR_ITERATION && delta_p_corrector>PRES_TOL);
    // delete node3
    DeleteNode(node3);
    //Save the solution properties 
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"P",n4.P);
    SetNodeData(node4,"rho",n4.rho);
    SetNodeData(node4,"T",n4.T);
    SetNodeData(node4,"V",n4.V);
    SetNodeData(node4,"theta",n4.theta);
    //connect the node into the characteristic mesh
    if (n4.x>n1.x){
	SetNodeData(node4,"CMinusUp",node1);
	SetNodeData(node1,"CMinusDown",node4);
    } else {
	SetNodeData(node4,"CMinusDown",node1);
	SetNodeData(node1,"CMinusUp",node4);
    } // end if
    if (n4.x>n2.x){
	SetNodeData(node4,"CPlusUp",node2);
	SetNodeData(node2,"CPlusDown",node4);
    } else {
	SetNodeData(node4,"CPlusDown",node2);
	SetNodeData(node2,"CPlusUp",node4);
    } // end if
    //Assuming a successful calculation, return the index of the solution node
    return node4;
} // end InteriorNode()
unittest {
    SetIsentropicFlag(NO); // homentropic solution
    SetFlowState(0); // work in 0th flow state
    SetAxiFlag(YES); // axisymmetric solution
    SetRatioOfSpecificHeats(1.20);
    SetGasConstant(320.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    //
    CreateNode(0);
    SetNodeData(0,"X",0.131460);
    SetNodeData(0,"Y",0.040118);
    SetNodeData(0,"V",2603.5);
    SetNodeData(0,"theta",18.191*PI/180.);
    SetNodeData(0,"P",34042.0);
    SetNodeData(0,"rho",0.086151);
    SetNodeData(0,"T",34042.0/R/0.086151);
    //
    CreateNode(1);
    SetNodeData(1,"X",0.135683);
    SetNodeData(1,"Y",0.037123);
    SetNodeData(1,"V",2609.2);
    SetNodeData(1,"theta",16.422*PI/180.);
    SetNodeData(1,"P",32781.0);
    SetNodeData(1,"rho",0.083482);
    SetNodeData(1,"T",32781.0/R/0.083482);
    //
    int sol = InteriorNode(0,1,-1);
    auto n4 = GetNodeData(sol);
    assert(approxEqual(n4.x,0.14118), "InteriorNode unittest: x-position failure");
    assert(approxEqual(n4.y,0.040554), "InteriorNode unittest: y-position failure");
    assert(approxEqual(n4.P,28767.), "InteriorNode unittest: pressure failure");
    assert(approxEqual(n4.rho,0.074875), "InteriorNode unittest: density failure");
    assert(approxEqual(n4.T,1200.63), "InteriorNode unittest: temperature failure");
    assert(approxEqual(n4.V,2628.6), "InteriorNode unittest: velocity failure");
    assert(approxEqual(n4.theta,17.267*PI/180.), "InteriorNode unittest: theta failure");
    DeleteNode(0);DeleteNode(1);DeleteNode(sol);
} // end InteriorNode unittest

/**
 * Purpose: Calculate a wall point from one initial (C+) point
 * and another point already on the wall
 * Input:
 * iw: Index of selected wall 
 *   node0: index of initial point on wall
 *   node2: index of initial point along node0's C- characteristic 
 *   node4: index of solution point (may have a value of -1) 
 * If -1 is specifed as the index for node4, a new node will be
 * created for the solution point.
 *   fs: flow state to work in
 * Output:
 * Returns the index of the solution point 
 */
int CPlusWallNode(int iw,int node0,int node2,int node4,int fs=-1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    //
    if(node0==-1 && GetIsentropicFlag==YES){
	return CPlusWallNode_0(iw,node2,node4,fs);
    } // end if
    if(node0 == node2){
	throw new Error(text("CPlusWallNode: node0 & node2 have same index: ",node2));
    } // end if
    if(CheckSameLine(node0,node2,-1) != YES){
	throw new Error(text("CPlusWallNode: node2 is not on C- char. of node0"));
    } // end if
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //first, assign the Node data
    NodeData n0 = GetNodeData(node0);
    if (n0 is null){
	throw new Error(text("CPlusWallNode: node0 doesn't exist, id: ",node0));
    } // end if
    if(CheckPointOnWall(iw,n0.pos) != YES){
	throw new Error(text("CPlusWallNode: node0, id:",node0," is not on wall, iw: ",iw));
    } // end if
    NodeData n2 = GetNodeData(node2);
    if (n2 is null){
	throw new Error(text("CPlusWallNode: node2 doesn't exist, id: ",node2));
    } // end if
    if (node4 == -1){
	node4 = CreateNode(-1);
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("StepStreamNode: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    // define gas properties
    double mu2 = MachAngle(n2.Mach);
    double mu0 = MachAngle(n0.Mach);
    //predictor variables
    double lambdaPlus = tan(n2.theta+mu2);
    double QPlus = sqrt(n2.Mach^^2-1.)/n2.rho/(n2.V^^2);
    double SPlus = axiflag*sin(n2.theta)/n2.y/n2.Mach/cos(n2.theta+mu2);
    //location of point 4, for predictor
    double t = WallFindT(iw,n2.pos,cos(atan(lambdaPlus)),sin(atan(lambdaPlus)));
    n4.pos = WallPos(iw,t);
    //
    double TPlus = -SPlus*(n4.x-n2.x)+QPlus*n2.P+n2.theta;
    double lambda0 = WallSlope(iw,t);
    n4.theta = atan(lambda0);
    //flow properties of point 4, for predictor
    double R0 = n0.rho*n0.V;
    double A0 = n0.a^^2;
    double Tee01 = R0*n0.V + n0.P;
    double Tee02 = n0.P - A0*n0.rho;
    //
    n4.P = (TPlus-n4.theta)/QPlus;
    n4.V = (Tee01-n4.P)/R0;
    n4.rho = (n4.P-Tee02)/A0;
    n4.T = n4.P/n4.rho/R;
    //solve for corrector
    int count=0;
    double P4_old,delta_p_corrector=0.0;
    do{
	count++;
	P4_old = n4.P;
	//coefficients, for corrector
	double pPlus = 0.5*(n2.P+n4.P);
	double thetaPlus = 0.5*(n2.theta+n4.theta);
	double VPlus = 0.5*(n2.V+n4.V);
	double rhoPlus = 0.5*(n2.rho+n4.rho);
	double yPlus = 0.5*(n2.y+n4.y);
	double aPlus = sqrt(g*pPlus/rhoPlus);
	double MPlus = VPlus/aPlus;
	double muPlus = MachAngle(MPlus);
	lambdaPlus = tan(thetaPlus+muPlus);
	QPlus = sqrt(MPlus^^2-1.)/rhoPlus/VPlus^^2;
	SPlus = axiflag*sin(thetaPlus)/yPlus/MPlus/cos(thetaPlus+muPlus);
	//location of point 4, for corrector
	t = WallFindT(iw,n2.pos,cos(atan(lambdaPlus)),sin(atan(lambdaPlus)));
	n4.pos = WallPos(iw,t);
	//
	TPlus = -SPlus*(n4.x-n2.x)+QPlus*n2.P+n2.theta;
	lambda0 = WallSlope(iw,t);
	n4.theta = atan(lambda0);
	//flow properties of point 4, for corrector
	double p0 = 0.5*(n0.P+n4.P);
	double rho0 = 0.5*(n0.rho+n4.rho);
	double a0 = sqrt(g*p0/rho0);
	double V0 = 0.5*(n0.V+n4.V);
	R0 = rho0*V0;
	A0 = a0^^2;
	Tee01 = R0*V0 + p0;
	Tee02 = p0 - A0*rho0;
	//
	n4.P = (TPlus-n4.theta)/QPlus;
	n4.V = (Tee01-n4.P)/R0;
	n4.rho = (n4.P-Tee02)/A0;
	n4.T = n4.P/n4.rho/R;
	//
	delta_p_corrector = fabs(n4.P - P4_old);
    } while (count<CORRECTOR_ITERATION && delta_p_corrector>PRES_TOL);
    //Save the solution properties 
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"P",n4.P);
    SetNodeData(node4,"rho",n4.rho);
    SetNodeData(node4,"T",n4.T);
    SetNodeData(node4,"V",n4.V);
    SetNodeData(node4,"theta",n4.theta);
    //connect the node into the characteristic mesh
    if (n4.x>n2.x){
	SetNodeData(node4,"CPlusUp",node2);
	SetNodeData(node2,"CPlusDown",node4);
    } else {
	SetNodeData(node4,"CPlusDown",node2);
	SetNodeData(node2,"CPlusUp",node4);
    } // end if
    //Assuming a successful calculation, return the index of the solution node
    return node4;
} // end CPlusWallNode()
unittest {

    SetIsentropicFlag(NO); // isentropic solution
    SetFlowState(0); // work in 0th flow state
    SetAxiFlag(YES); // axisymmetric solution
    SetRatioOfSpecificHeats(1.20);
    SetGasConstant(320.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    //
    WallFromPolynomial(0,[0.0221852,0.71568,-1.0787],0.05,0.08);
    //
    CreateNode(0);
    SetNodeData(0,"X",0.055943);
    SetNodeData(0,"Y",0.058845);
    SetNodeData(0,"V",2252.9);
    SetNodeData(0,"theta",30.752*PI/180.);
    SetNodeData(0,"P",2.1453e5);
    SetNodeData(0,"rho",0.39947);
    SetNodeData(0,"T",2.1453e5/R/0.39947);
    SetNodeData(0,"CMD",1);
    //
    CreateNode(1);
    SetNodeData(1,"X",0.060480);
    SetNodeData(1,"Y",0.059625);
    SetNodeData(1,"V",2274.2);
    SetNodeData(1,"theta",30.122*PI/180.);
    SetNodeData(1,"P",1.9597e5);
    SetNodeData(1,"rho",0.37046);
    SetNodeData(1,"T",1.9597e5/R/0.37046);
    SetNodeData(1,"CMU",0);
    //
    int sol = CPlusWallNode(0,0,1,-1);
    auto n4 = GetNodeData(sol);
    //
    assert(approxEqual(n4.x,0.063486), "CPlusWallNode unittest: x-position failure");
    assert(approxEqual(n4.y,0.063272), "CPlusWallNode unittest: y-position failure");
    assert(approxEqual(n4.P,1.8729e5), "CPlusWallNode unittest: pressure failure");
    assert(approxEqual(n4.rho,0.35675), "CPlusWallNode unittest: density failure");
    assert(approxEqual(n4.T,1.8729e5/R/0.35675), "CPlusWallNode unittest: temperature failure");
    assert(approxEqual(n4.V,2284.7), "CPlusWallNode unittest: velocity failure");
    assert(approxEqual(n4.theta,30.058*PI/180.), "CPlusWallNode unittest: theta failure");
    DeleteWall(0);DeleteNode(0);DeleteNode(1);DeleteNode(sol);
} // end CPlusWallNode unittest

/**
 * Purpose: Calculate a wall point from one initial (C-) point
 * and another point already on the wall
 * Input:
 * iw: Index of selected wall 
 *   node0: index of initial point on wall
 *   node1: index of initial point along node0's C+ characteristic 
 *   node4: index of solution point (may have a value of -1) 
 * If -1 is specifed as the index for node4, a new node will be
 * created for the solution point.
 *   fs: flow state to work in
 * Output:
 * Returns the index of the solution point 
 */
int CMinusWallNode(int iw,int node0,int node1,int node4,int fs=-1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    //
    if(node0==-1 && GetIsentropicFlag==YES){
	return CMinusWallNode_0(iw,node1,node4,fs);
    } // end if
    if(node0 == node1){
	throw new Error(text("CMinusWallNode: node0 & node1 have same index: ",node1));
    } // end if
    if(CheckSameLine(node0,node1,1) != YES){
	throw new Error(text("CMinusWallNode: node1 is not on C+ char. of node0"));
    } // end if
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //first, assign the Node data
    NodeData n0 = GetNodeData(node0);
    if (n0 is null){
	throw new Error(text("CMinusWallNode: node0 doesn't exist, id: ",node0));
    } // end if
    if(CheckPointOnWall(iw,n0.pos) != YES){
	throw new Error(text("CMinusWallNode: node0, id:",node0," is not on wall, iw: ",iw));
    } // end if
    NodeData n1 = GetNodeData(node1);
    if (n1 is null){
	throw new Error(text("CMinusWallNode: node1 doesn't exist, id: ",node1));
    } // end if
    if (node4 == -1){
	node4 = CreateNode(-1);
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("StepStreamNode: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    // define gas properties
    double mu1 = MachAngle(n1.Mach);
    double mu0 = MachAngle(n0.Mach);
    //predictor variables
    double lambdaMinus = tan(n1.theta-mu1);
    double QMinus = sqrt(n1.Mach^^2-1.)/n1.rho/(n1.V^^2);
    double SMinus = axiflag*sin(n1.theta)/n1.y/n1.Mach/cos(n1.theta-mu1);
    //location of point 4, for predictor
    double t = WallFindT(iw,n1.pos,cos(atan(lambdaMinus)),sin(atan(lambdaMinus)));
    n4.pos = WallPos(iw,t);
    //
    double TMinus = -SMinus*(n4.x-n1.x)+QMinus*n1.P-n1.theta;
    double lambda0 = WallSlope(iw,t);
    n4.theta = atan(lambda0);
    //flow properties of point 4, for predictor
    double R0 = n0.rho*n0.V;
    double A0 = n0.a^^2;
    double Tee01 = R0*n0.V + n0.P;
    double Tee02 = n0.P - A0*n0.rho;
    //
    n4.P = (TMinus+n4.theta)/QMinus;
    n4.V = (Tee01-n4.P)/R0;
    n4.rho = (n4.P-Tee02)/A0;
    n4.T = n4.P/n4.rho/R;
    //solve for corrector
    int count=0;
    double P4_old,delta_p_corrector=0.0;
    do{
	count++;
	P4_old = n4.P;
	//coefficients, for corrector
	double pMinus = 0.5*(n1.P+n4.P);
	double thetaMinus = 0.5*(n1.theta+n4.theta);
	double VMinus = 0.5*(n1.V+n4.V);
	double rhoMinus = 0.5*(n1.rho+n4.rho);
	double yMinus = 0.5*(n1.y+n4.y);
	double aMinus = sqrt(g*pMinus/rhoMinus);
	double MMinus = VMinus/aMinus;
	double muMinus = MachAngle(MMinus);
	lambdaMinus = tan(thetaMinus-muMinus);
	QMinus = sqrt(MMinus^^2-1.)/rhoMinus/VMinus^^2;
	SMinus = axiflag*sin(thetaMinus)/yMinus/MMinus/cos(thetaMinus-muMinus);
	//location of point 4, for corrector
	t = WallFindT(iw,n1.pos,cos(atan(lambdaMinus)),sin(atan(lambdaMinus)));
	n4.pos = WallPos(iw,t);
	//
	TMinus = -SMinus*(n4.x-n1.x)+QMinus*n1.P-n1.theta;
	lambda0 = WallSlope(iw,t);
	n4.theta = atan(lambda0);
	//flow properties of point 4, for corrector
	double p0 = 0.5*(n0.P+n4.P);
	double rho0 = 0.5*(n0.rho+n4.rho);
	double a0 = sqrt(g*p0/rho0);
	double V0 = 0.5*(n0.V+n4.V);
	R0 = rho0*V0;
	A0 = a0^^2;
	Tee01 = R0*V0 + p0;
	Tee02 = p0 - A0*rho0;
	//
	n4.P = (TMinus+n4.theta)/QMinus;
	n4.V = (Tee01-n4.P)/R0;
	n4.rho = (n4.P-Tee02)/A0;
	n4.T = n4.P/n4.rho/R;
	//
	delta_p_corrector = fabs(n4.P - P4_old);
    } while (count<CORRECTOR_ITERATION && delta_p_corrector>PRES_TOL);
    //Save the solution properties 
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"P",n4.P);
    SetNodeData(node4,"rho",n4.rho);
    SetNodeData(node4,"T",n4.T);
    SetNodeData(node4,"V",n4.V);
    SetNodeData(node4,"theta",n4.theta);
    //connect the node into the characteristic mesh
    if (n4.x>n1.x){
	SetNodeData(node4,"CMinusUp",node1);
	SetNodeData(node1,"CMinusDown",node4);
    } else {
	SetNodeData(node4,"CMinusDown",node1);
	SetNodeData(node1,"CMinusUp",node4);
    } // end if
    //Assuming a successful calculation, return the index of the solution node
    return node4;
} // end CMinusWallNode()
unittest {
    SetIsentropicFlag(NO); // isentropic solution
    SetFlowState(0); // work in 0th flow state
    SetAxiFlag(YES); // axisymmetric solution
    SetRatioOfSpecificHeats(1.20);
    SetGasConstant(320.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    //
    WallFromPolynomial(0,[-0.0221852,-0.71568,1.0787],0.05,0.08);
    //
    CreateNode(0);
    SetNodeData(0,"X",0.055943);
    SetNodeData(0,"Y",-0.058845);
    SetNodeData(0,"V",2252.9);
    SetNodeData(0,"theta",-30.752*PI/180.);
    SetNodeData(0,"P",2.1453e5);
    SetNodeData(0,"rho",0.39947);
    SetNodeData(0,"T",2.1453e5/R/0.39947);
    SetNodeData(0,"CPD",1);
    //
    CreateNode(1);
    SetNodeData(1,"X",0.060480);
    SetNodeData(1,"Y",-0.059625);
    SetNodeData(1,"V",2274.2);
    SetNodeData(1,"theta",-30.122*PI/180.);
    SetNodeData(1,"P",1.9597e5);
    SetNodeData(1,"rho",0.37046);
    SetNodeData(1,"T",1.9597e5/R/0.37046);
    SetNodeData(1,"CPU",0);
    //
    int sol = CMinusWallNode(0,0,1,-1);
    auto n4 = GetNodeData(sol);
    //
    assert(approxEqual(n4.x,0.063486), "CMinusWallNode unittest: x-position failure");
    assert(approxEqual(n4.y,-0.063272), "CMinusWallNode unittest: y-position failure");
    assert(approxEqual(n4.P,1.8729e5), "CMinusWallNode unittest: pressure failure");
    assert(approxEqual(n4.rho,0.35675), "CMinusWallNode unittest: density failure");
    assert(approxEqual(n4.T,1.8729e5/R/0.35675), "CMinusWallNode unittest: temperature failure");
    assert(approxEqual(n4.V,2284.7), "CMinusWallNode unittest: velocity failure");
    assert(approxEqual(n4.theta,-30.058*PI/180.), "CMinusWallNode unittest: theta failure");
    DeleteWall(0);DeleteNode(0);DeleteNode(1);DeleteNode(sol);
} // end CMinusWallNode unittest

/**
 * Purpose: Calculate a free-boundary point from 
 * one point (node0) already on the boundary and 
 * one point (node2) on a C+ characteristic.
 * Input: 
 *   node0: index of initial point along C0 streamline
 *   node2: index of initial point along C+ characteristic
 *   node4: index of solution point (may have a value of -1)
 * If -1 is specified as the index for node4, a new node will be
 * created for the solution point.
 *   fs: flow state to work in
 * Output:
 * Returns the index of the solution point
 */
int CPlusFreeBndyNode(int node0,int node2,int node4,int fs=-1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    //
    if(GetIsentropicFlag==YES){
	return CPlusFreeBndyNode_0(node0,node2,node4,fs);
    } // end if
    if(node0 == node2){
	throw new Error(text("CPlusFreeBndyNode: node0 & node2 have same index: ",node2));
    } // end if
    if(CheckSameLine(node0,node2,-1) != YES){
	throw new Error(text("CPlusFreeBndyNode: node2 is not on C- char. of node0"));
    } // end if
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //first, assign the Node data
    NodeData n0 = GetNodeData(node0);
    if (n0 is null){
	throw new Error(text("CPlusFreeBndyNode: node0 doesn't exist, id: ",node0));
    } // end if
    NodeData n2 = GetNodeData(node2);
    if (n2 is null){
	throw new Error(text("CPlusFreeBndyNode: node2 doesn't exist, id: ",node2));
    } // end if
    if (node4 == -1){
	node4 = CreateNode(-1);
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("StepStreamNode: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    //gas properties along free boundary are constant
    n4.P = n0.P;
    n4.rho = n0.rho;
    n4.T = n0.T;
    n4.V = n0.V;
    // define gas properties
    double mu2 = MachAngle(n2.Mach);
    double mu0 = MachAngle(n0.Mach);
    //predictor variables
    double lambdaPlus = tan(n2.theta+mu2);
    double QPlus = sqrt(n2.Mach^^2-1.)/n2.rho/(n2.V^^2);
    double SPlus = axiflag*sin(n2.theta)/n2.y/n2.Mach/cos(n2.theta+mu2);
    double lambda0 = tan(n0.theta);
    //location of point 4, for predictor
    n4.x = (n0.y-n2.y-lambdaPlus*n2.x+lambda0*n0.x)/(lambda0-lambdaPlus);
    n4.y = n0.y-lambda0*(n0.x-n4.x);
    //
    double TPlus = -SPlus*(n4.x-n2.x)+QPlus*n2.P+n2.theta;
    //
    n4.theta = TPlus - QPlus*n4.P;
    //solve for corrector
    double pPlus = 0.5*(n2.P+n4.P);
    double VPlus = 0.5*(n2.V+n4.V);
    double rhoPlus = 0.5*(n2.rho+n4.rho);
    //
    Vector3 xy4;
    int count=0;
    double delta_position;
    do{
	count++;
	xy4 = n4.pos;
	//coefficients, for corrector
	double thetaPlus = 0.5*(n2.theta+n4.theta);
	//
	double yPlus = 0.5*(n2.y+n4.y);
	double aPlus = sqrt(g*pPlus/rhoPlus);
	double MPlus = VPlus/aPlus;
	double muPlus = MachAngle(MPlus);
	lambdaPlus = tan(thetaPlus+muPlus);
	QPlus = sqrt(MPlus^^2-1.)/rhoPlus/VPlus^^2;
	SPlus = axiflag*sin(thetaPlus)/yPlus/MPlus/cos(thetaPlus+muPlus);
	//location of point 4, for corrector
	n4.x = (n0.y-n2.y-lambdaPlus*n2.x+lambda0*n0.x)/(lambda0-lambdaPlus);
	n4.y = n0.y-lambda0*(n0.x-n4.x);
	//
	TPlus = -SPlus*(n4.x-n2.x)+QPlus*n2.P+n2.theta;
	//
	n4.theta = TPlus - QPlus*n4.P;
	//
	delta_position = abs(xy4-n4.pos);
    }while (count<CORRECTOR_ITERATION && delta_position>POS_TOL);
    //Save the solution properties 
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"P",n4.P);
    SetNodeData(node4,"rho",n4.rho);
    SetNodeData(node4,"T",n4.T);
    SetNodeData(node4,"V",n4.V);
    SetNodeData(node4,"theta",n4.theta);
    //connect the node into the characteristic mesh
    if (n4.x>n2.x){
	SetNodeData(node4,"CPlusUp",node2);
	SetNodeData(node2,"CPlusDown",node4);
    } else {
	SetNodeData(node4,"CPlusDown",node2);
	SetNodeData(node2,"CPlusUp",node4);
    } // end if
    //Assuming a successful calculation, return the index of the solution node
    return node4;
} // end CPlusFreeBndyNode()
unittest {
    SetIsentropicFlag(NO); // isentropic solution
    SetFlowState(0); // work in 0th flow state
    SetAxiFlag(YES); // axisymmetric solution
    SetRatioOfSpecificHeats(1.20);
    SetGasConstant(320.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    //
    CreateNode(0);
    SetNodeData(0,"X",0.32979);
    SetNodeData(0,"Y",0.12351);
    SetNodeData(0,"V",2511.7);
    SetNodeData(0,"theta",15.317*PI/180.);
    SetNodeData(0,"P",60000.);
    SetNodeData(0,"T",60000./R/0.13813);
    SetNodeData(0,"rho",0.13813);
    SetNodeData(0,"CMD",1);
    //
    CreateNode(1);
    SetNodeData(1,"X",0.34226);
    SetNodeData(1,"Y",0.12312);
    SetNodeData(1,"V",2532.30);
    SetNodeData(1,"theta",14.165*PI/180.);
    SetNodeData(1,"P",53164.0);
    SetNodeData(1,"T",53164.0/R/0.12491);
    SetNodeData(1,"rho",0.12491);
    SetNodeData(1,"CMU",0);
    //
    int sol = CPlusFreeBndyNode(0,1,-1);
    auto n4 = GetNodeData(sol);
    //
    assert(approxEqual(n4.x,0.35282), "CPlusFreeBndyNode unittest: x-position failure");
    assert(approxEqual(n4.y,0.12916), "CPlusFreeBndyNode unittest: y-position failure");
    assert(approxEqual(n4.P,6.0e4), "CPlusFreeBndyNode unittest: pressure failure");
    assert(approxEqual(n4.rho,0.13813), "CPlusFreeBndyNode unittest: density failure");
    assert(approxEqual(n4.T,60000./R/0.13813), "CPlusFreeBndyNode unittest: temperature failure");
    assert(approxEqual(n4.V,2511.7), "CPlusFreeBndyNode unittest: velocity failure");
    assert(approxEqual(n4.theta,12.230*PI/180.), "CPlusFreeBndyNode unittest: theta failure");
    DeleteNode(0);DeleteNode(1);DeleteNode(2);DeleteNode(sol);
} // end CPlusFreeBndyNode unittest

/**
 * Purpose: Calculate a free-boundary point from 
 * one point (node0) already on the boundary and 
 * one point (node1) on a C- characteristic.
 * Input: 
 *   node0: index of initial point along C0 streamline
 *   node1: index of initial point along C- characteristic
 *   node4: index of solution point (may have a value of -1)
 * If -1 is specified as the index for node4, a new node will be
 * created for the solution point.
 *   fs: flow state to work in
 * Output:
 * Returns the index of the solution point
 */
int CMinusFreeBndyNode(int node0,int node1,int node4,int fs=-1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    //
    if(GetIsentropicFlag==YES){
	return CMinusFreeBndyNode_0(node1,node4,fs);
    } // end if
    if(node0 == node1){
	throw new Error(text("CMinusFreeBndyNode: node0 & node1 have same index: ",node1));
    } // end if
    if(CheckSameLine(node0,node1,1) != YES){
	throw new Error(text("CPlusFreeBndyNode: node1 is not on C+ char. of node0"));
    } // end if
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //first, assign the Node data
    NodeData n0 = GetNodeData(node0);
    if (n0 is null){
	throw new Error(text("CMinusFreeBndyNode: node0 doesn't exist, id: ",node0));
    } // end if
    NodeData n1 = GetNodeData(node1);
    if (n1 is null){
	throw new Error(text("CMinusFreeBndyNode: node1 doesn't exist, id: ",node1));
    } // end if
    if (node4 == -1){
	node4 = CreateNode(-1);
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("StepStreamNode: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    //gas properties along free boundary are constant
    n4.P = n0.P;
    n4.rho = n0.rho;
    n4.T = n0.T;
    n4.V = n0.V;
    // define gas properties
    double mu1 = MachAngle(n1.Mach);
    double mu0 = MachAngle(n0.Mach);
    //predictor variables
    double lambdaMinus = tan(n1.theta-mu1);
    double QMinus = sqrt(n1.Mach^^2-1.)/n1.rho/(n1.V^^2);
    double SMinus = axiflag*sin(n1.theta)/n1.y/n1.Mach/cos(n1.theta-mu1);
    double lambda0 = tan(n0.theta);
    //location of point 4, for predictor
    n4.x = (n0.y-n1.y-lambdaMinus*n1.x+lambda0*n0.x)/(lambda0-lambdaMinus);
    n4.y = n0.y-lambda0*(n0.x-n4.x);
    //
    double TMinus = -SMinus*(n4.x-n1.x)+QMinus*n1.P-n1.theta;
    //
    n4.theta = -TMinus + QMinus*n4.P;
    //solve for corrector
    double pMinus = 0.5*(n1.P+n4.P);
    double VMinus = 0.5*(n1.V+n4.V);
    double rhoMinus = 0.5*(n1.rho+n4.rho);
    //
    Vector3 xy4;
    int count=0;
    double delta_position;
    do{
	count++;
	xy4 = n4.pos;
	//coefficients, for corrector
	double thetaMinus = 0.5*(n1.theta+n4.theta);
	//
	double yMinus = 0.5*(n1.y+n4.y);
	double aMinus = sqrt(g*pMinus/rhoMinus);
	double MMinus = VMinus/aMinus;
	double muMinus = MachAngle(MMinus);
	lambdaMinus = tan(thetaMinus-muMinus);
	QMinus = sqrt(MMinus^^2-1.)/rhoMinus/VMinus^^2;
	SMinus = axiflag*sin(thetaMinus)/yMinus/MMinus/cos(thetaMinus-muMinus);
	lambda0 = tan(0.5*(n0.theta+n4.theta));
	//location of point 4, for corrector
	n4.x = (n0.y-n1.y-lambdaMinus*n1.x+lambda0*n0.x)/(lambda0-lambdaMinus);
	n4.y = n0.y-lambda0*(n0.x-n4.x);
	//
	TMinus = -SMinus*(n4.x-n1.x)+QMinus*n1.P-n1.theta;
	//
	n4.theta = -TMinus + QMinus*n4.P;
	//
	delta_position = abs(xy4-n4.pos);
    }while (count<CORRECTOR_ITERATION && delta_position>POS_TOL);
    //Save the solution properties 
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"P",n4.P);
    SetNodeData(node4,"rho",n4.rho);
    SetNodeData(node4,"T",n4.T);
    SetNodeData(node4,"V",n4.V);
    SetNodeData(node4,"theta",n4.theta);
    //connect the node into the characteristic mesh
    if (n4.x>n1.x){
	SetNodeData(node4,"CMinusUp",node1);
	SetNodeData(node1,"CMinusDown",node4);
    } else {
	SetNodeData(node4,"CMinusDown",node1);
	SetNodeData(node1,"CMinusUp",node4);
    } // end if
    //Assuming a successful calculation, return the index of the solution node
    return node4;
} // end CMinusFreeBndyNode()
unittest {

    SetIsentropicFlag(YES); // isentropic solution
    SetFlowState(0); // work in 0th flow state
    SetAxiFlag(YES); // axisymmetric solution
    SetRatioOfSpecificHeats(1.20);
    SetGasConstant(320.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    //
    CreateNode(0);
    SetNodeData(0,"X",0.32979);
    SetNodeData(0,"Y",-0.12351);
    SetNodeData(0,"V",2511.72);
    SetNodeData(0,"theta",-15.317*PI/180.);
    SetNodeData(0,"P",60000.);
    SetNodeData(0,"rho",0.13813);
    SetNodeData(0,"T",60000./R/0.13813);
    SetNodeData(0,"CPD",1);
    //
    CreateNode(1);
    SetNodeData(1,"X",0.34226);
    SetNodeData(1,"Y",-0.12312);
    SetNodeData(1,"V",2532.30);
    SetNodeData(1,"theta",-14.165*PI/180.);
    SetNodeData(1,"P",53164.);
    SetNodeData(1,"rho",0.12491);
    SetNodeData(1,"T",53164./R/0.12491);
    SetNodeData(1,"CPU",0);
    //
    int sol = CMinusFreeBndyNode_0(0,1,-1);
    auto n4 = GetNodeData(sol);
    //
    assert(approxEqual(n4.x,0.35283), "CMinusFreeBndyNode unittest: x-position failure");
    assert(approxEqual(n4.y,-0.12916), "CMinusFreeBndyNode unittest: y-position failure");
    assert(approxEqual(n4.P,6.0e4), "CMinusFreeBndyNode unittest: pressure failure");
    assert(approxEqual(n4.rho,6.0e4/R/1357.14), "CMinusFreeBndyNode unittest: density failure");
    assert(approxEqual(n4.T,1357.14), "CMinusFreeBndyNode unittest: temperature failure");
    assert(approxEqual(n4.V,2511.7), "CMinusFreeBndyNode unittest: velocity failure");
    assert(approxEqual(n4.theta,-12.23*PI/180.), "CMinusFreeBndyNode unittest: theta failure");
    DeleteNode(0);DeleteNode(1);DeleteNode(sol);
} // end CMinusFreeBndyNode unittest

/**
 * Purpose: Calculate a new post-shock point on an oblique shock wave from 
 * one initial shock point (node0) and a point on its C- char. (node2).
 * The solution node will be on the C+ char. of node2.
 * Input:
 *   iw: Index of selected wall. 
 *   node0: index of initial point oblique shock
 *   node2: index of initial point on node0's C- line
 *   node4: index of solution point (may have a value of -1) 
 *   If -1 is specified as the index for node4, a new node will be
 *   created for the solution point.
 *   fs: flow state to work in
 * Output :
 *   Returns the index of the solution point
 */

int CPlusOblShkNode(int node0,int node2,int node4,int fs=1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    fs = GetFlowState();
	
    if (CheckSameLine(node0,node2,-1) != 1){
	throw new Error(text("CPlusOblShkNode: node0 and node2 are not on the same C- char. line"));
    } // end if

    if(node0 == node2){
	throw new Error(text("CPlusOblShkNode: node0 & node2 have same index: ",node2));
    } // end if
    //first, assign the Node data
    NodeData n0 = GetNodeData(node0);
    if (n0 is null){
	throw new Error(text("CPlusOblShkNode: node0 doesn't exist, id: ",node0));
    } // end if
    if (GetNumberFlowStates(node0) < 2){
	throw new Error(text("CPlusOblShkNode: node0 only has 1 flowstate"));
    } // end if
    NodeData n2 = GetNodeData(node2);
    if (n2 is null){
	throw new Error(text("CPlusOblShkNode: node2 doesn't exist, id: ",node2));
    } // end if
    if (node4 == -1){
	node4 = CreateNode(-1);
    }
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("StepStreamNode: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    n4.F[fs-1].T = n0.F[fs-1].T;
    n4.F[fs-1].P = n0.F[fs-1].P;
    n4.F[fs-1].rho = n0.F[fs-1].rho;
    n4.F[fs-1].V = n0.F[fs-1].V;
    n4.F[fs-1].theta = n0.F[fs-1].theta;
    MakeOblShkNode(node4,"theta",n4.F[fs-1].theta,fs);
    //
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //initial guess at solution point, will be very wrong but used to check convergence
    double beta4 = n0.beta(fs);
    //compute actual solution point
    int iteration_count=0;
    double dP1,dP2;
    double beta_old1,beta_old2;
    double delta_pressure;
    do {
	++iteration_count;
	if (iteration_count>=3){
	    beta4 = beta_old2+(beta_old2-beta_old1)*dP2/(dP1-dP2);
	    beta_old1 = beta_old2; dP1 = dP2;
	    beta_old2 = beta4; dP2 = delta_pressure;
	} // end if
	double lambda0 = tan(0.5*(n0.beta+beta4));
	n4.theta = n4.F[fs-1].theta + theta_obl(n4.F[fs-1].Mach,beta4,g);
	n4.P = n4.F[fs-1].P*p2_p1_obl(n4.F[fs-1].Mach,beta4,g);
	n4.V = n4.F[fs-1].V*V2_V1_obl(n4.F[fs-1].Mach,beta4,n4.theta-n4.F[fs-1].theta,g);
	n4.rho = n4.F[fs-1].rho*r2_r1_obl(n4.F[fs-1].Mach,beta4,g);
	n4.T = n4.P/R/n4.rho;
	//
	double thetaPlus = 0.5*(n2.theta + n4.theta); 
	double MPlus = 0.5*(n2.Mach + n4.Mach);
	double pPlus = 0.5*(n2.P+n4.P);
	double muPlus = MachAngle(MPlus);
	double lambdaPlus = tan(muPlus+thetaPlus);
	double x4 = (lambdaPlus*n2.x-lambda0*n0.x+n0.y-n2.y)/(lambdaPlus-lambda0);
	double y4 = lambdaPlus*(x4-n2.x) + n2.y;
	n4.pos = Vector3(x4,y4);
	double QPlus = sqrt(MPlus^^2-1.0)/g/pPlus/(MPlus^^2);
	double SPlus = axiflag*sin(thetaPlus)/y4/MPlus/cos(thetaPlus+muPlus);
	double TPlus = -SPlus*(x4-n2.x) + QPlus*n2.P + n2.theta;
	double P4_LRC = (TPlus-n4.theta)/QPlus;

	delta_pressure = P4_LRC - n4.P;
	if (iteration_count==1){
	    dP1 = delta_pressure;
	    beta4 = beta_obl2(n4.F[fs-1].Mach,P4_LRC/n4.F[fs-1].P,g);
	    beta_old1 = beta4;
	} // end if
	if (iteration_count==2){
	    dP2 = delta_pressure;
	    beta4 = beta_obl2(n4.F[fs-1].Mach,P4_LRC/n4.F[fs-1].P,g);
	    beta_old2 = beta4;
	} // end if
	if (iteration_count>=3){
	    dP2 = delta_pressure;
	} // end if
    } while(iteration_count<MAX_ITERATION && fabs(delta_pressure)>PRES_TOL);
    //Save the solution-point properties
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"P",n4.P);
    SetNodeData(node4,"rho",n4.rho);
    SetNodeData(node4,"T",n4.T);
    SetNodeData(node4,"V",n4.V);
    SetNodeData(node4,"theta",n4.theta);
    //connect the node into the characteristic mesh.
    if (n4.x>n2.x){
	SetNodeData(node4,"CPlusUp",node2);
	SetNodeData(node2,"CPlusDown",node4);
    } else {
	SetNodeData(node4,"CPlusDown",node2);
	SetNodeData(node2,"CPlusUp",node4);
    } // end if
    return node4;
} // end CPlusOblShkNode

unittest {

    SetFlowState(1); // work in 0th flow state
    SetAxiFlag(NO); // axisymmetric solution
    SetRatioOfSpecificHeats(1.40);
    SetGasConstant(287.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    //
    CreateNode(0);
    SetNodeData(0,"X",0.22113,0);
    SetNodeData(0,"Y",0.28317,0);
    SetNodeData(0,"T",300.0,0);
    SetNodeData(0,"P",1.0e5,0);
    SetNodeData(0,"rho",1.0e5/R/300.0,0);
    SetNodeData(0,"theta",0.0,0);
    SetNodeData(0,"V", sqrt(g*R*300.) * 3.0,0);
    SetNodeData(0,"CMD",1,0);
    //
    MakeOblShkNode(0,"theta",30.*PI/180.,1);
    //
    CreateNode(1);
    SetNodeData(1,"X",0.26991);
    SetNodeData(1,"Y",0.26988);
    SetNodeData(1,"V",697.16);
    SetNodeData(1,"theta",29.536*PI/180.);
    SetNodeData(1,"P",6.2132e5);
    SetNodeData(1,"rho",3.6189);
    SetNodeData(1,"T",6.2132e5/R/3.6189);
    SetNodeData(1,"CMU",0);
    //
    int sol = CPlusOblShkNode(0,1,-1,1);
    auto n4 = GetNodeData(sol);
    //
    assert(approxEqual(n4.x,0.30376), "CPlusOblShkNode unittest: x-position failure");
    assert(approxEqual(n4.y,0.38744), "CPlusOblShkNode unittest: y-position failure");
    assert(approxEqual(n4.P,6.2095e5), "CPlusOblShkNode unittest: pressure failure");
    assert(approxEqual(n4.rho,3.6387), "CPlusOblShkNode unittest: density failure");
    assert(approxEqual(n4.T,6.2095e5/R/3.6387), "CPlusOblShkNode unittest: temperature failure");
    assert(approxEqual(n4.V,702.31), "CPlusOblShkNode unittest: velocity failure");
    assert(approxEqual(n4.theta,29.548*PI/180.), "CPlusOblShkNode unittest: theta failure");
    DeleteNode(0);DeleteNode(1);DeleteNode(sol);
} // end CPlusOblShkNode unittest

/**
 * Purpose: Calculate a new post-shock point on an oblique shock wave from 
 * one initial shock point (node0) and a point on its C+ char. (node2).
 * The solution node will be on the C- char. of node1.
 * Input:
 *   iw: Index of selected wall. 
 *   node0: index of initial point oblique shock
 *   node1: index of initial point on node0's C+ line
 *   node4: index of solution point (may have a value of -1) 
 *   If -1 is specified as the index for node4, a new node will be
 *   created for the solution point.
 *   fs: flow state to work in
 * Output :
 *   Returns the index of the solution point
 */

int CMinusOblShkNode(int node0,int node1,int node4,int fs=1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    fs = GetFlowState();
    if (CheckSameLine(node0,node1,1) != 1){
	throw new Error(text("CMinusOblShkNode: node0 and node1 are not on the same C- char. line"));
    } // end if
    if(node0 == node1){
	throw new Error(text("CMinusOblShkNode: node0 & node1 have same index: ",node1));
    } // end if
    //first, assign the Node data
    NodeData n0 = GetNodeData(node0);
    if (n0 is null){
	throw new Error(text("CMinusOblShkNode: node0 doesn't exist, id: ",node0));
    } // end if
    if (GetNumberFlowStates(node0) < 2){
	throw new Error(text("CMinusOblShkNode: node0 only has 1 flowstate"));
    } // end if
    NodeData n1 = GetNodeData(node1);
    if (n1 is null){
	throw new Error(text("CMinusOblShkNode: node1 doesn't exist, id: ",node1));
    } // end if
    if (node4 == -1){
	node4 = CreateNode(-1);
    }
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("StepStreamNode: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    n4.F[fs-1].T = n0.F[fs-1].T;
    n4.F[fs-1].P = n0.F[fs-1].P;
    n4.F[fs-1].rho = n0.F[fs-1].rho;
    n4.F[fs-1].V = n0.F[fs-1].V;
    n4.F[fs-1].theta = n0.F[fs-1].theta;
    MakeOblShkNode(node4,"theta",n4.F[fs-1].theta,fs);
    //
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //initial guess at solution point, will be very wrong but used to check convergence
    double beta4 = n0.beta();
    //compute actual solution point
    int iteration_count=0;
    double dP1,dP2;
    double beta_old1,beta_old2;
    double delta_pressure;
    do {
	++iteration_count;
	if (iteration_count>=3){
	    beta4 = beta_old2+(beta_old2-beta_old1)*dP2/(dP1-dP2);
	    beta_old1 = beta_old2; dP1 = dP2;
	    beta_old2 = beta4; dP2 = delta_pressure;
	} // end if
	double lambda0 = tan(0.5*(n0.beta+beta4));
	n4.theta = n4.F[fs-1].theta + theta_obl(n4.F[fs-1].Mach,beta4,g);
	n4.P = n4.F[fs-1].P*p2_p1_obl(n4.F[fs-1].Mach,beta4,g);
	n4.V = n4.F[fs-1].V*V2_V1_obl(n4.F[fs-1].Mach,beta4,n4.theta-n4.F[fs-1].theta,g);
	n4.rho = n4.F[fs-1].rho*r2_r1_obl(n4.F[fs-1].Mach,beta4,g);
	n4.T = n4.P/R/n4.rho;
	//
	double thetaMinus = 0.5*(n1.theta + n4.theta); 
	double MMinus = 0.5*(n1.Mach + n4.Mach);
	double pMinus = 0.5*(n1.P+n4.P);
	double muMinus = MachAngle(MMinus);
	double lambdaMinus = tan(thetaMinus-muMinus);
	double x4 = (lambdaMinus*n1.x-lambda0*n0.x+n0.y-n1.y)/(lambdaMinus-lambda0);
	double y4 = lambdaMinus*(x4-n1.x) + n1.y;
	n4.pos = Vector3(x4,y4);
	double QMinus = sqrt(MMinus^^2-1.0)/g/pMinus/(MMinus^^2);
	double SMinus = axiflag*sin(thetaMinus)/y4/MMinus/cos(thetaMinus-muMinus);
	double TMinus = -SMinus*(x4-n1.x) + QMinus*n1.P - n1.theta;
	double P4_LRC = (TMinus+n4.theta)/QMinus;

	delta_pressure = P4_LRC - n4.P;
	if (iteration_count==1){
	    dP1 = delta_pressure;
	    beta4 = -beta_obl2(n4.F[fs-1].Mach,P4_LRC/n4.F[fs-1].P,g);
	    beta_old1 = beta4;
	} // end if
	if (iteration_count==2){
	    dP2 = delta_pressure;
	    beta4 = -beta_obl2(n4.F[fs-1].Mach,P4_LRC/n4.F[fs-1].P,g);
	    beta_old2 = beta4;
	} // end if
	if (iteration_count>=3){
	    dP2 = delta_pressure;
	} // end if
    } while(iteration_count<MAX_ITERATION && fabs(delta_pressure)>PRES_TOL);
    //Save the solution-point properties
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"P",n4.P);
    SetNodeData(node4,"rho",n4.rho);
    SetNodeData(node4,"T",n4.T);
    SetNodeData(node4,"V",n4.V);
    SetNodeData(node4,"theta",n4.theta);
    //connect the node into the characteristic mesh.
    if (n4.x>n1.x){
	SetNodeData(node4,"CMinusUp",node1);
	SetNodeData(node1,"CMinusDown",node4);
    } else{
	SetNodeData(node4,"CMinusDown",node1);
	SetNodeData(node1,"CMinusUp",node4);
    } // end if
    return node4;
} // end CMinusOblShkNode()
unittest {
    SetFlowState(1); // work in 0th flow state
    SetAxiFlag(NO); // axisymmetric solution
    SetRatioOfSpecificHeats(1.40);
    SetGasConstant(287.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    //
    CreateNode(0);
    SetNodeData(0,"X",0.22113,0);
    SetNodeData(0,"Y",-0.28317,0);
    SetNodeData(0,"T",300.0,0);
    SetNodeData(0,"P",1.0e5,0);
    SetNodeData(0,"rho",1.0e5/R/300.0,0);
    SetNodeData(0,"theta",0.0,0);
    SetNodeData(0,"V", sqrt(g*R*300.) * 3.0,0);
    SetNodeData(0,"CPD",1,0);
    //
    MakeOblShkNode(0,"theta",-30.*PI/180.,1);
    //
    CreateNode(1);
    SetNodeData(1,"X",0.26991);
    SetNodeData(1,"Y",-0.26988);
    SetNodeData(1,"V",697.16);
    SetNodeData(1,"theta",-29.536*PI/180.);
    SetNodeData(1,"P",6.2132e5);
    SetNodeData(1,"rho",3.6189);
    SetNodeData(1,"T",6.2132e5/R/3.6189);
    SetNodeData(1,"CPU",0,0);
    //
    int sol = CMinusOblShkNode(0,1,-1,1);
    auto n4 = GetNodeData(sol);
    //
    assert(approxEqual(n4.x,0.30376), "CPlusOblShkNode unittest: x-position failure");
    assert(approxEqual(n4.y,-0.38744), "CPlusOblShkNode unittest: y-position failure");
    assert(approxEqual(n4.P,6.2095e5), "CPlusOblShkNode unittest: pressure failure");
    assert(approxEqual(n4.rho,3.6387), "CPlusOblShkNode unittest: density failure");
    assert(approxEqual(n4.T,6.2095e5/R/3.6387), "CPlusOblShkNode unittest: temperature failure");
    assert(approxEqual(n4.V,702.31), "CPlusOblShkNode unittest: velocity failure");
    assert(approxEqual(n4.theta,-29.548*PI/180.), "CPlusOblShkNode unittest: theta failure");
    DeleteNode(0);DeleteNode(1);DeleteNode(sol);	
} // end CPlusOblShkNode unittest

//************************************************************//
// Isentropic unit processes                                  //
// all from Zucrow M, Hoffman J, 1985, Gas Dynamics:          //
// Multidimensional Flow, Volume II, John Wiley and Sons      //
//************************************************************//

/**
 * Calculate an interior point from two initial points
 * Input:
 *   node1: index of initial point along C- characteristic
 *   node2: index of initial point along C+ characteristic
 *   node4: index of solution point (may have a value of -1)
 * If -1 is specified as the index for node4, a new node will be
 * created for the solution point.
 *   fs: flow state to work in
 * Output: 
 * Returns the index of the solution point 
 */
int InteriorNode_0(int node1,int node2,int node4,int fs=-1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    //
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //first, assign the Node data
    if(node1 == node2){
	throw new Error(text("InteriorNode_0: node1 & node2 have same id: ",node1));
    } // end if
    NodeData n1 = GetNodeData(node1);
    if (n1 is null){
	throw new Error(text("InteriorNode_0: node1 doesn't exist, id: ",node1));
    } // end if
    NodeData n2 = GetNodeData(node2);
    if (n2 is null){
	throw new Error(text("InteriorNode_0: node2 doesn't exist, id: ",node2));
    } // end if
    if (node4 == -1){
	node4 = CreateNode(-1);
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("InteriorNode_0: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    // define gas properties
    double mu1 = MachAngle(n1.Mach);
    double mu2 = MachAngle(n2.Mach);
    double u1 = n1.V*cos(n1.theta), v1 = n1.V*sin(n1.theta);
    double u2 = n2.V*cos(n2.theta), v2 = n2.V*sin(n2.theta);
    // predictor variables
    double lambdaPlus = tan(n2.theta+mu2);
    double QPlus = u2^^2-n2.a^^2;
    double RPlus = 2.*u2*v2 - QPlus*lambdaPlus;
    double SPlus = axiflag*n2.a^^2*v2/n2.y;
    //
    double lambdaMinus = tan(n1.theta-mu1);
    double QMinus = u1^^2-n1.a^^2;
    double RMinus = 2.*u1*v1 - QMinus*lambdaMinus;
    double SMinus = axiflag*n1.a^^2*v1/n1.y;
    //location of point 4, for predictor
    n4.x = (n1.y-n2.y-lambdaMinus*n1.x+lambdaPlus*n2.x)/(lambdaPlus-lambdaMinus);
    n4.y = n2.y-lambdaPlus*(n2.x-n4.x);
    //
    double TPlus = SPlus*(n4.x-n2.x)+QPlus*u2+RPlus*v2;
    double TMinus = SMinus*(n4.x-n1.x)+QMinus*u1+RMinus*v1;
    double v4 = (TMinus-TPlus*QMinus/QPlus)/(RMinus-RPlus*QMinus/QPlus);
    double u4 = (TPlus-RPlus*v4)/QPlus;
    //
    double x4_old,y4_old;
    int iteration_count=0;
    double delta_position=0.0;
    do {
	++iteration_count;
	x4_old = n4.x;y4_old = n4.y;
	//
	double uPlus = 0.5*(u2+u4);
	double vPlus = 0.5*(v2+v4);
	double yPlus = 0.5*(n2.y+n4.y);
	//
	double uMinus = 0.5*(u1+u4);
	double vMinus = 0.5*(v1+v4);
	double yMinus = 0.5*(n1.y+n4.y);
	//
	double VPlus = sqrt(uPlus^^2+vPlus^^2);
	double thetaPlus = atan(vPlus /uPlus);
	double aPlus = sqrt(g*R*n2.T0 - 0.5*(g-1.)*VPlus^^2);
	double muPlus = asin(aPlus/VPlus);
	//
	double VMinus = sqrt(uMinus^^2+vMinus^^2);
	double thetaMinus = atan(vMinus /uMinus);
	double aMinus = sqrt(g*R*n1.T0 - 0.5*(g-1.)*VMinus^^2);
	double muMinus = asin(aMinus/VMinus);
	//
	lambdaPlus = tan(thetaPlus+muPlus);
	QPlus = uPlus^^2-aPlus^^2;
	RPlus = 2.*uPlus*vPlus - QPlus*lambdaPlus;
	SPlus = axiflag*aPlus^^2*vPlus/yPlus;
	//
	lambdaMinus = tan(thetaMinus-muMinus);
	QMinus = uMinus^^2-aMinus^^2;
	RMinus = 2.*uMinus*vMinus - QMinus*lambdaMinus;
	SMinus = axiflag*aMinus^^2*vMinus/yMinus;
	//location of point 4, for corrector
	n4.x = (n1.y-n2.y-lambdaMinus*n1.x+lambdaPlus*n2.x)/(lambdaPlus-lambdaMinus);
	n4.y = n2.y-lambdaPlus*(n2.x-n4.x);
	//
	TPlus = SPlus*(n4.x-n2.x)+QPlus*u2+RPlus*v2;
	TMinus = SMinus*(n4.x-n1.x)+QMinus*u1+RMinus*v1;
	v4 = (TMinus-TPlus*QMinus/QPlus)/(RMinus-RPlus*QMinus/QPlus);
	u4 = (TPlus-RPlus*v4)/QPlus;
	//
	delta_position=sqrt((n4.x-x4_old)^^2+(n4.y-y4_old)^^2);                 
    } while (iteration_count<MAX_ITERATION && delta_position>POS_TOL);
    //Save the solution properties 
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"V",sqrt(u4^^2+v4^^2));
    SetNodeData(node4,"T",n1.T0 - n4.V^^2*(g-1.)/2./g/R);
    SetNodeData(node4,"P",n1.P0 / p0_p(n4.Mach,g));
    SetNodeData(node4,"rho",n4.P/R/n4.T);
    SetNodeData(node4,"theta",atan(v4/u4));
    //connect the node into the characteristic mesh
    if (n4.x>n1.x){
	SetNodeData(node4,"CMinusUp",node1);
	SetNodeData(node1,"CMinusDown",node4);
    } else {
	SetNodeData(node4,"CMinusDown",node1);
	SetNodeData(node1,"CMinusUp",node4);
    } // end if
    if (n4.x>n2.x){
	SetNodeData(node4,"CPlusUp",node2);
	SetNodeData(node2,"CPlusDown",node4);
    } else {
	SetNodeData(node4,"CPlusDown",node2);
	SetNodeData(node2,"CPlusUp",node4);
    } // end if
    //Assuming a successful calculation, return the index of the solution node
    return node4;
} // InteriorNode_0()
unittest {
    SetIsentropicFlag(YES); // isentropic solution
    SetFlowState(0); // work in 0th flow state
    SetAxiFlag(YES); // axisymmetric solution
    SetRatioOfSpecificHeats(1.20);
    SetGasConstant(320.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    //
    CreateNode(0);
    SetNodeData(0,"X",0.131460);
    SetNodeData(0,"Y",0.040118);
    SetNodeData(0,"V",2603.5);
    SetNodeData(0,"theta",18.191*PI/180.);
    SetNodeData(0,"P",34042.0);
    SetNodeData(0,"rho",0.086151);
    SetNodeData(0,"T",34042.0/R/0.086151);
    //
    CreateNode(1);
    SetNodeData(1,"X",0.135683);
    SetNodeData(1,"Y",0.037123);
    SetNodeData(1,"V",2609.2);
    SetNodeData(1,"theta",16.422*PI/180.);
    SetNodeData(1,"P",32781.0);
    SetNodeData(1,"rho",0.083482);
    SetNodeData(1,"T",32781.0/R/0.083482);
    //
    int sol = InteriorNode_0(0,1,-1);
    auto n4 = GetNodeData(sol);
    //
    assert(approxEqual(n4.x,0.14118), "InteriorNode_0 unittest: x-position failure");
    assert(approxEqual(n4.y,0.040554), "InteriorNode_0 unittest: y-position failure");
    assert(approxEqual(n4.P,28767.), "InteriorNode_0 unittest: pressure failure");
    assert(approxEqual(n4.rho,0.074875), "InteriorNode_0 unittest: density failure");
    assert(approxEqual(n4.T,1200.63), "InteriorNode_0 unittest: temperature failure");
    assert(approxEqual(n4.V,2628.6), "InteriorNode_0 unittest: velocity failure");
    assert(approxEqual(n4.theta,17.267*PI/180.), "InteriorNode_0 unittest: theta failure");
    DeleteNode(0);DeleteNode(1);DeleteNode(sol);
} // end InteriorNode_0 unittest

/**
 * Purpose: Calculate a wall point from one initial (C+) point.
 * Input:
 *   iw: Index of selected wall. 
 *   node2: index of initial point along C+ characteristic 
 *   node4: index of solution point (may have a value of -1) 
 * If -1 is specified as the index for node4, a new node will be
 * created for the solution point.
 *   fs: flow state to work in
 * Output:
 *   Returns the index of the solution point
 */
int CPlusWallNode_0(int iw,int node2,int node4,int fs=-1)
{
    // check wall exists
    if(WallIsPresent(iw) != YES){
	throw new Error(text("CPlusWallNode_0: Wall %d is not present.",iw));
    } // end if
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    //
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //first, assign the Node data
    NodeData n2 = GetNodeData(node2);
    if (n2 is null){
	throw new Error(text("CPlusWallNode_0: node2 doesn't exist, id: ",node2));
    } // end if
    if (node4 == -1){
	node4 = CreateNode(-1);
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("StepStreamNode: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    // define gas properties
    double mu2 = MachAngle(n2.Mach);
    double u2 = n2.V*cos(n2.theta);
    double v2 = n2.V*sin(n2.theta);
    // predictor variables
    double lambdaPlus = tan(n2.theta+mu2);
    double QPlus = u2^^2-n2.a^^2;
    double RPlus = 2.*u2*v2 - QPlus*lambdaPlus;
    double SPlus = axiflag*n2.a^^2*v2/n2.y;
    //location of point 4, for predictor
    double t = WallFindT(iw,n2.pos,cos(atan(lambdaPlus)),sin(atan(lambdaPlus)));
    n4.pos = WallPos(iw,t);
    //
    double TPlus = SPlus*(n4.x-n2.x)+QPlus*u2+RPlus*v2;
    double lambda0 = WallSlope(iw,t);
    n4.theta = atan(lambda0);
    double u4 = TPlus/(QPlus+lambda0*RPlus);
    double v4 = u4*lambda0;
    //
    Vector3 pos4_old;
    int iteration_count=0;
    double delta_position=0.0;
    do {
	++iteration_count;
	pos4_old = n4.pos;
	//
	double uPlus = 0.5*(u2+u4);
	double vPlus = 0.5*(v2+v4);
	double yPlus = 0.5*(n2.y+n4.y);
	//
	double VPlus = sqrt(uPlus^^2+vPlus^^2);
	double thetaPlus = atan(vPlus /uPlus);
	double aPlus = sqrt(g*R*n2.T0 - 0.5*(g-1.)*VPlus^^2);
	double muPlus = asin(aPlus/VPlus);
	//
	lambdaPlus = tan(thetaPlus+muPlus);
	QPlus = uPlus^^2-aPlus^^2;
	RPlus = 2.*uPlus*vPlus - QPlus*lambdaPlus;
	SPlus = axiflag*aPlus^^2*vPlus/yPlus;
	//location of point 4, for corrector
	t = WallFindT(iw,n2.pos,cos(atan(lambdaPlus)),sin(atan(lambdaPlus)));
	n4.pos = WallPos(iw,t);
	//
	TPlus = SPlus*(n4.x-n2.x)+QPlus*u2+RPlus*v2;
	lambda0 = WallSlope(iw,t);
	n4.theta = atan(lambda0);
	u4 = TPlus/(QPlus+lambda0*RPlus);
	v4 = u4*lambda0;
	//
	delta_position = abs(pos4_old - n4.pos);              
    } while (iteration_count<MAX_ITERATION && delta_position>POS_TOL);
    //Save the solution properties 
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"V",sqrt(u4^^2+v4^^2));
    SetNodeData(node4,"T",n2.T0 - n4.V^^2*(g-1.)/2./g/R);
    SetNodeData(node4,"P",n2.P0 / p0_p(n4.Mach,g));
    SetNodeData(node4,"rho",n4.P/R/n4.T);
    SetNodeData(node4,"theta",n4.theta);
    //connect the node into the characteristic mesh
    if (n4.x>n2.x){
	SetNodeData(node4,"CPlusUp",node2);
	SetNodeData(node2,"CPlusDown",node4);
    }
    else{
	SetNodeData(node4,"CPlusDown",node2);
	SetNodeData(node2,"CPlusUp",node4);
    } // end if
    //Assuming a successful calculation, return the index of the solution node.
    return node4;
} // end CPlusWallNode_0()

unittest {

    SetIsentropicFlag(YES); // isentropic solution
    SetFlowState(0); // work in 0th flow state
    SetAxiFlag(YES); // axisymmetric solution
    SetRatioOfSpecificHeats(1.20);
    SetGasConstant(320.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    //
    WallFromPolynomial(0,[0.0221852,0.71568,-1.0787],0.05,0.08);
    //
    CreateNode(0);
    SetNodeData(0,"X",0.060480);
    SetNodeData(0,"Y",0.059625);
    SetNodeData(0,"V",2274.2);
    SetNodeData(0,"theta",30.12*PI/180.);
    SetNodeData(0,"P",2.37440e5);
    SetNodeData(0,"rho",2.37440e5/R/1653.13);
    SetNodeData(0,"T",1653.13);
    //
    int sol = CPlusWallNode_0(0,0,-1);
    auto n4 = GetNodeData(sol);
    //
    assert(approxEqual(n4.x,0.063485), "CPlusWallNode_0 unittest: x-position failure");
    assert(approxEqual(n4.y,0.063273), "CPlusWallNode_0 unittest: y-position failure");
    assert(approxEqual(n4.P,2.29071e5), "CPlusWallNode_0 unittest: pressure failure");
    assert(approxEqual(n4.rho,2.29071e5/R/1640.7), "CPlusWallNode_0 unittest: density failure");
    assert(approxEqual(n4.T,1640.7), "CPlusWallNode_0 unittest: temperature failure");
    assert(approxEqual(n4.V,2284.7), "CPlusWallNode_0 unittest: velocity failure");
    assert(approxEqual(n4.theta,30.06*PI/180.), "CPlusWallNode_0 unittest: theta failure");
    DeleteWall(0);DeleteNode(0);DeleteNode(1);
} // end CPlusWallNode_0 unittest

/**
 * Purpose: Calculate a wall point from one initial (C-) point.
 * Input:
 *   iw: Index of selected wall. 
 *   node1: index of initial point along C- characteristic 
 *   node4: index of solution point (may have a value of -1) 
 * If -1 is specified as the index for node4, a new node will be
 * created for the solution point.
 *   fs: flow state to work in
 * Output:
 *   Returns the index of the solution point
 */
int CMinusWallNode_0(int iw,int node1,int node4,int fs=-1)
{
    // check wall exists
    if(WallIsPresent(iw) != YES){
	throw new Error(text("CMinusWallNode_0: Wall %d is not present.",iw));
    } // end if
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    //
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //first, assign the Node data
    NodeData n1 = GetNodeData(node1);
    if (n1 is null){
	throw new Error(text("CMinusWallNode_0: node1 doesn't exist, id: ",node1));
    } // end if
    if (node4 == -1){
	node4 = CreateNode(-1);
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("StepStreamNode: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    // define gas properties
    double mu1 = MachAngle(n1.Mach);
    double u1 = n1.V*cos(n1.theta);
    double v1 = n1.V*sin(n1.theta);
    // predictor variables
    double lambdaMinus = tan(n1.theta-mu1);
    double QMinus = u1^^2-n1.a^^2;
    double RMinus = 2.*u1*v1 - QMinus*lambdaMinus;
    double SMinus = axiflag*n1.a^^2*v1/n1.y;
    //location of point 4, for predictor
    double t = WallFindT(iw,n1.pos,cos(atan(lambdaMinus)),sin(atan(lambdaMinus)));
    n4.pos = WallPos(iw,t);
    //
    double TMinus = SMinus*(n4.x-n1.x)+QMinus*u1+RMinus*v1;
    double lambda0 = WallSlope(iw,t);
    n4.theta = atan(lambda0);
    double u4 = TMinus/(QMinus+lambda0*RMinus);
    double v4 = u4*lambda0;
    //
    Vector3 pos4_old;
    int iteration_count=0;
    double delta_position=0.0;
    do {
	++iteration_count;
	pos4_old = n4.pos;
	//
	double uMinus = 0.5*(u1+u4);
	double vMinus = 0.5*(v1+v4);
	double yMinus = 0.5*(n1.y+n4.y);
	//
	double VMinus = sqrt(uMinus^^2+vMinus^^2);
	double thetaMinus = atan(vMinus/uMinus);
	double aMinus = sqrt(g*R*n1.T0 - 0.5*(g-1.)*VMinus^^2);
	double muMinus = asin(aMinus/VMinus);
	//
	lambdaMinus = tan(thetaMinus-muMinus);
	QMinus = uMinus^^2-aMinus^^2;
	RMinus = 2.*uMinus*vMinus - QMinus*lambdaMinus;
	SMinus = axiflag*aMinus^^2*vMinus/yMinus;
	//location of point 4, for corrector
	t = WallFindT(iw,n1.pos,cos(atan(lambdaMinus)),sin(atan(lambdaMinus)));
	n4.pos = WallPos(iw,t);
	//
	TMinus = SMinus*(n4.x-n1.x)+QMinus*u1+RMinus*v1;
	lambda0 = WallSlope(iw,t);
	n4.theta = atan(lambda0);
	u4 = TMinus/(QMinus+lambda0*RMinus);
	v4 = u4*lambda0;
	//
	delta_position = abs(pos4_old - n4.pos);              
    } while (iteration_count<MAX_ITERATION && delta_position>POS_TOL);
    //Save the solution properties 
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"V",sqrt(u4^^2+v4^^2));
    SetNodeData(node4,"T",n1.T0 - n4.V^^2*(g-1.)/2./g/R);
    SetNodeData(node4,"P",n1.P0 / p0_p(n4.Mach,g));
    SetNodeData(node4,"rho",n4.P/R/n4.T);
    SetNodeData(node4,"theta",n4.theta);
    //connect the node into the characteristic mesh
    if (n4.x>n1.x){
	SetNodeData(node4,"CMinusUp",node1);
	SetNodeData(node1,"CMinusDown",node4);
    }
    else{
	SetNodeData(node4,"CMinusDown",node1);
	SetNodeData(node1,"CMinusUp",node4);
    } // end if
    //Assuming a successful calculation, return the index of the solution node.
    return node4;
} // end CMinusWallNode_0()

unittest {

    SetIsentropicFlag(YES); // isentropic solution
    SetFlowState(0); // work in 0th flow state
    SetAxiFlag(YES); // axisymmetric solution
    SetRatioOfSpecificHeats(1.20);
    SetGasConstant(320.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    //
    WallFromPolynomial(0,[-0.0221852,-0.71568,1.0787],0.05,0.08);
    //
    CreateNode(0);
    SetNodeData(0,"X",0.060480);
    SetNodeData(0,"Y",-0.059625);
    SetNodeData(0,"V",2274.2);
    SetNodeData(0,"theta",-30.12*PI/180.);
    SetNodeData(0,"P",2.37440e5);
    SetNodeData(0,"rho",2.37440e5/R/1653.13);
    SetNodeData(0,"T",1653.13);
    //
    int sol = CMinusWallNode_0(0,0,-1);
    auto n4 = GetNodeData(sol);
    //
    assert(approxEqual(n4.x,0.063485), "CMinusWallNode_0 unittest: x-position failure");
    assert(approxEqual(n4.y,-0.063273), "CMinusWallNode_0 unittest: y-position failure");
    assert(approxEqual(n4.P,2.29071e5), "CMinusWallNode_0 unittest: pressure failure");
    assert(approxEqual(n4.rho,2.29071e5/R/1640.7), "CMinusWallNode_0 unittest: density failure");
    assert(approxEqual(n4.T,1640.7), "CMinusWallNode_0 unittest: temperature failure");
    assert(approxEqual(n4.V,2284.7), "CMinusWallNode_0 unittest: velocity failure");
    assert(approxEqual(n4.theta,-30.06*PI/180.), "CMinusWallNode_0 unittest: theta failure");
    DeleteWall(0);DeleteNode(0);DeleteNode(1);
} // end CMinusWallNode_0 unittest

/**
 * Purpose: Calculate a free-boundary point from 
 * one point (node0) already on the boundary and 
 * one point (node2) on a C+ characteristic.
 * Input: 
 *   node0: index of initial point along C0 streamline
 *   node2: index of initial point along C+ characteristic
 *   node4: index of solution point (may have a value of -1)
 * If -1 is specified as the index for node4, a new node will be
 * created for the solution point.
 *   fs: flow state to work in
 * Output:
 * Returns the index of the solution point
 */
int CPlusFreeBndyNode_0(int node0,int node2,int node4,int fs=-1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    if(node2 == node0){
	throw new Error(text("CPlusFreeBndyNode: node0 & node2 have same index: ",node2));
    } // end if
    if(CheckSameLine(node0,node2,-1) != YES){
	throw new Error(text("CPlusFreeBndyNode: node2 is not on C- char. of node0"));
    } // end if
    //
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //first, assign the Node data
    NodeData n0 = GetNodeData(node0);
    if (n0 is null){
	throw new Error(text("CPlusFreeBndyNode_0: node0 doesn't exist, id: ",node0));
    } // end if
    NodeData n2 = GetNodeData(node2);
    if (n2 is null){
	throw new Error(text("CPlusFreeBndyNode_0: node2 doesn't exist, id: ",node2));
    } // end if
    if (node4 == -1){
	node4 = CreateNode(-1);
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("StepStreamNode: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    n4 = n0.dup;
    // define gas properties
    double mu2 = MachAngle(n2.Mach);
    double u2 = n2.V*cos(n2.theta);
    double v2 = n2.V*sin(n2.theta);
    int sign=1;
    if(n2.theta < 0.0){sign=-1;}
    // predictor variables
    double lambdaPlus = tan(n2.theta+mu2);
    double QPlus = u2^^2-n2.a^^2;
    double RPlus = 2.*u2*v2 - QPlus*lambdaPlus;
    double SPlus = axiflag*n2.a^^2*v2/n2.y;
    double lambda0 = tan(n0.theta);
    //location of point 4, for predictor
    n4.x = (n0.y-n2.y-lambdaPlus*n2.x+lambda0*n0.x)/(lambda0-lambdaPlus);
    n4.y = n0.y-lambda0*(n0.x-n4.x);
    //
    double TPlus = SPlus*(n4.x-n2.x)+QPlus*u2+RPlus*v2;
    double u4 = (QPlus*TPlus-RPlus*sqrt(n4.V^^2*(QPlus^^2+RPlus^^2)-TPlus^^2))/(QPlus^^2+RPlus^^2);
    double v4 = sqrt(n4.V^^2 - u4^^2) * sign;
    //
    double x4_old,y4_old;
    int iteration_count=0;
    double delta_position=0.0;
    do {
	++iteration_count;
	x4_old = n4.x;y4_old = n4.y;
	//
	double uPlus = 0.5*(u2+u4);
	double vPlus = 0.5*(v2+v4);
	double yPlus = 0.5*(n2.y+n4.y);
	//
	double VPlus = sqrt(uPlus^^2+vPlus^^2);
	double thetaPlus = atan(vPlus /uPlus);
	double aPlus = sqrt(g*R*n2.T0 - 0.5*(g-1.)*VPlus^^2);
	double muPlus = asin(aPlus/VPlus);
	//
	//
	lambdaPlus = tan(thetaPlus+muPlus);
	QPlus = uPlus^^2-aPlus^^2;
	RPlus = 2.*uPlus*vPlus - QPlus*lambdaPlus;
	SPlus = axiflag*aPlus^^2*vPlus/yPlus;
	lambda0 = 0.5*(tan(n0.theta)+v4/u4);
	//location of point 4, for corrector
	n4.x = (n0.y-n2.y-lambdaPlus*n2.x+lambda0*n0.x)/(lambda0-lambdaPlus);
	n4.y = n0.y-lambda0*(n0.x-n4.x);
	//
	TPlus = SPlus*(n4.x-n2.x)+QPlus*u2+RPlus*v2;
	u4 = (QPlus*TPlus-RPlus*sqrt(n4.V^^2*(QPlus^^2+RPlus^^2)-TPlus^^2))/(QPlus^^2+RPlus^^2);
	v4 = sqrt(n4.V^^2 - u4^^2) * sign;
	//
	delta_position=sqrt((n4.x-x4_old)^^2+(n4.y-y4_old)^^2);          
    } while (iteration_count<MAX_ITERATION && delta_position>POS_TOL);
    //Save the solution properties 
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"V",n4.V);
    SetNodeData(node4,"T",n4.T);
    SetNodeData(node4,"P",n4.P);
    SetNodeData(node4,"rho",n4.rho);
    SetNodeData(node4,"theta",atan(v4/u4));
    //connect the node into the characteristic mesh
    if (n4.x>n2.x){
	SetNodeData(node4,"CPlusUp",node2);
	SetNodeData(node2,"CPlusDown",node4);
    }
    else{
	SetNodeData(node4,"CPlusDown",node2);
	SetNodeData(node2,"CPlusUp",node4);
    } // end if
    //Assuming a successful calculation, return the index of the solution node.
    return node4;
} // end CPlusFreeBndyNode_0()

unittest {

    SetIsentropicFlag(YES); // isentropic solution
    SetFlowState(0); // work in 0th flow state
    SetAxiFlag(YES); // axisymmetric solution
    SetRatioOfSpecificHeats(1.20);
    SetGasConstant(320.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    //
    CreateNode(0);
    SetNodeData(0,"X",0.32979);
    SetNodeData(0,"Y",0.12351);
    SetNodeData(0,"V",2511.72);
    SetNodeData(0,"theta",15.317*PI/180.);
    SetNodeData(0,"P",60000.);
    SetNodeData(0,"rho",60000./R/1357.10);
    SetNodeData(0,"T",1357.10);
    SetNodeData(0,"CMD",1);
    //
    CreateNode(1);
    SetNodeData(1,"X",0.34226);
    SetNodeData(1,"Y",0.12312);
    SetNodeData(1,"V",2532.30);
    SetNodeData(1,"theta",14.165*PI/180.);
    SetNodeData(1,"P",53174.3);
    SetNodeData(1,"rho",53174.3/R/1330.07);
    SetNodeData(1,"T",1330.07);
    SetNodeData(1,"CMU",0);
    //
    int sol = CPlusFreeBndyNode_0(0,1,-1);
    auto n4 = GetNodeData(sol);
    //
    assert(approxEqual(n4.x,0.35283), "CPlusFreeBndyNode_0 unittest: x-position failure");
    assert(approxEqual(n4.y,0.12916), "CPlusFreeBndyNode_0 unittest: y-position failure");
    assert(approxEqual(n4.P,6.0e4), "CPlusFreeBndyNode_0 unittest: pressure failure");
    assert(approxEqual(n4.rho,6.0e4/R/1357.14), "CPlusFreeBndyNode_0 unittest: density failure");
    assert(approxEqual(n4.T,1357.14), "CPlusFreeBndyNode_0 unittest: temperature failure");
    assert(approxEqual(n4.V,2511.7), "CPlusFreeBndyNode_0 unittest: velocity failure");
    assert(approxEqual(n4.theta,12.23*PI/180.), "CPlusFreeBndyNode_0 unittest: theta failure");
    DeleteNode(0);DeleteNode(1);DeleteNode(sol);
} // end CPlusFreeBndyNode_0 unittest

/**
 * Purpose: Calculate a free-boundary point from 
 * one point (node0) already on the boundary and 
 * one point (node1) on a C- characteristic.
 * Input: 
 *   node0: index of initial point along C0 streamline
 *   node1: index of initial point along C- characteristic
 *   node4: index of solution point (may have a value of -1)
 * If -1 is specified as the index for node4, a new node will be
 * created for the solution point.
 *   fs: flow state to work in
 * Output:
 * Returns the index of the solution point
 */
int CMinusFreeBndyNode_0(int node0,int node1,int node4,int fs=-1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    if(node1 == node0){
	throw new Error(text("CMinusFreeBndyNode: node0 & node1 have same index: ",node1));
    } // end if
    if(CheckSameLine(node0,node1,1) != YES){
	throw new Error(text("CMinusFreeBndyNode: node1 is not on C+ char. of node0"));
    } // end if
    //
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //first, assign the Node data
    NodeData n0 = GetNodeData(node0);
    if (n0 is null){
	throw new Error(text("CMinusFreeBndyNode_0: node0 doesn't exist, id: ",node0));
    } // end if
    NodeData n1 = GetNodeData(node1);
    if (n1 is null){
	throw new Error(text("CMinusFreeBndyNode_0: node1 doesn't exist, id: ",node1));
    } // end if
    if (node4 == -1){
	node4 = CreateNode(-1);
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("StepStreamNode: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    n4 = n0.dup;
    // define gas properties
    double mu1 = MachAngle(n1.Mach);
    double u1 = n1.V*cos(n1.theta);
    double v1 = n1.V*sin(n1.theta);
    int sign=1;
    if(n1.theta < 0.0){sign=-1;}
    // predictor variables
    double lambdaMinus = tan(n1.theta-mu1);
    double QMinus = u1^^2-n1.a^^2;
    double RMinus = 2.*u1*v1 - QMinus*lambdaMinus;
    double SMinus = axiflag*n1.a^^2*v1/n1.y;
    double lambda0 = tan(n0.theta);
    //location of point 4, for predictor
    n4.x = (n0.y-n1.y-lambdaMinus*n1.x+lambda0*n0.x)/(lambda0-lambdaMinus);
    n4.y = n0.y-lambda0*(n0.x-n4.x);
    //
    double TMinus = SMinus*(n4.x-n1.x)+QMinus*u1+RMinus*v1;
    double u4 = (QMinus*TMinus+RMinus*sqrt(n4.V^^2*(QMinus^^2+RMinus^^2)-TMinus^^2))/(QMinus^^2+RMinus^^2);
    double v4 = sqrt(n4.V^^2 - u4^^2) * sign;
    //
    double x4_old,y4_old;
    int iteration_count=0;
    double delta_position=0.0;
    do {
	++iteration_count;
	x4_old = n4.x;y4_old = n4.y;
	//
	double uMinus = 0.5*(u1+u4);
	double vMinus = 0.5*(v1+v4);
	double yMinus = 0.5*(n1.y+n4.y);
	//
	double VMinus = sqrt(uMinus^^2+vMinus^^2);
	double thetaMinus = atan(vMinus /uMinus);
	double aMinus = sqrt(g*R*n1.T0 - 0.5*(g-1.)*VMinus^^2);
	double muMinus = asin(aMinus/VMinus);
	//
	//
	lambdaMinus = tan(thetaMinus-muMinus);
	QMinus = uMinus^^2-aMinus^^2;
	RMinus = 2.*uMinus*vMinus - QMinus*lambdaMinus;
	SMinus = axiflag*aMinus^^2*vMinus/yMinus;
	lambda0 = 0.5*(tan(n0.theta)+v4/u4);
	//location of point 4, for corrector
	n4.x = (n0.y-n1.y-lambdaMinus*n1.x+lambda0*n0.x)/(lambda0-lambdaMinus);
	n4.y = n0.y-lambda0*(n0.x-n4.x);
	//
	TMinus = SMinus*(n4.x-n1.x)+QMinus*u1+RMinus*v1;
	u4 = (QMinus*TMinus+RMinus*sqrt(n4.V^^2*(QMinus^^2+RMinus^^2)-TMinus^^2))/(QMinus^^2+RMinus^^2);
	v4 = sqrt(n4.V^^2 - u4^^2) * sign;
	//
	delta_position=sqrt((n4.x-x4_old)^^2+(n4.y-y4_old)^^2);  
    } while (iteration_count<MAX_ITERATION && delta_position>POS_TOL);
    //Save the solution properties 
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"V",n4.V);
    SetNodeData(node4,"T",n4.T);
    SetNodeData(node4,"P",n4.P);
    SetNodeData(node4,"rho",n4.rho);
    SetNodeData(node4,"theta",atan(v4/u4));
    //connect the node into the characteristic mesh
    if (n4.x>n1.x){
	SetNodeData(node4,"CMinusUp",node1);
	SetNodeData(node1,"CMinusDown",node4);
    }
    else{
	SetNodeData(node4,"CMinusDown",node1);
	SetNodeData(node1,"CMinusUp",node4);
    } // end if
    //Assuming a successful calculation, return the index of the solution node.
    return node4;
} // end CMinusFreeBndyNode_0()

unittest {

    SetIsentropicFlag(YES); // isentropic solution
    SetFlowState(0); // work in 0th flow state
    SetAxiFlag(YES); // axisymmetric solution
    SetRatioOfSpecificHeats(1.20);
    SetGasConstant(320.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    //
    CreateNode(0);
    SetNodeData(0,"X",0.32979);
    SetNodeData(0,"Y",-0.12351);
    SetNodeData(0,"V",2511.72);
    SetNodeData(0,"theta",-15.317*PI/180.);
    SetNodeData(0,"P",60000.);
    SetNodeData(0,"rho",60000./R/1357.10);
    SetNodeData(0,"T",1357.10);
    SetNodeData(0,"CPD",1);
    //
    CreateNode(1);
    SetNodeData(1,"X",0.34226);
    SetNodeData(1,"Y",-0.12312);
    SetNodeData(1,"V",2532.30);
    SetNodeData(1,"theta",-14.165*PI/180.);
    SetNodeData(1,"P",53174.3);
    SetNodeData(1,"rho",53174.3/R/1330.07);
    SetNodeData(1,"T",1330.07);
    SetNodeData(1,"CPU",0);
    //
    int sol = CMinusFreeBndyNode_0(0,1,2);
    auto n4 = GetNodeData(sol);
    //
    assert(approxEqual(n4.x,0.35283), "CMinusFreeBndyNode_0 unittest: x-position failure");
    assert(approxEqual(n4.y,-0.12916), "CMinusFreeBndyNode_0 unittest: y-position failure");
    assert(approxEqual(n4.P,6.0e4), "CMinusFreeBndyNode_0 unittest: pressure failure");
    assert(approxEqual(n4.rho,6.0e4/R/1357.14), "CMinusFreeBndyNode_0 unittest: density failure");
    assert(approxEqual(n4.T,1357.14), "CMinusFreeBndyNode_0 unittest: temperature failure");
    assert(approxEqual(n4.V,2511.7), "CMinusFreeBndyNode_0 unittest: velocity failure");
    assert(approxEqual(n4.theta,-12.23*PI/180.), "CMinusFreeBndyNode_0 unittest: theta failure");
    DeleteNode(0);DeleteNode(1);DeleteNode(sol);
} // end CMinusFreeBndyNode_0 unittest

//************************************************************//
// Other processes                                            //
//************************************************************//

/**
 * Purpose: Generate nodes along a new C-/+ characteristic line
 * Input: 
 *   old_first: a starting point on an existing C-/+ line 
 *   new_first: the starting point on the new C-/+ line 
 *             on which the new nodes are to be generated
 *   char_line: plus or minus [Plus/plus/+ or Minus/minus/-]
 *   dir: up or down ["Up"/"up" or "Down"/"down"] 
 *   fs: flow state to work in
 * Output : 
 * Returns a list of node indices on the new C-/+ curve.
 * in the order that they are generated.
 */
int MarchAlongC(out int[] nodemarchlist,int old_first,int new_first,
		string char_line="minus",string dir="down",int fs=-1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    if(old_first == new_first){
	throw new Error(text("MarchAlongC: old_first & new_first have same index: ",old_first));
    } // end if
    NodeData n0 = GetNodeData(old_first);
    if (n0 is null){
	throw new Error(text("MarchAlongC: old_first node doesn't exist, id: ",old_first));
    } // end if
    NodeData n1 = GetNodeData(new_first);
    if (n1 is null){
	throw new Error(text("MarchAlongC: new_first node doesn't exist, id: ",new_first));
    } // end if
    //
    int char_flag;
    if(char_line=="minus" || char_line=="Minus" || char_line=="-"){
	char_flag = -1;
    } else if(char_line=="plus" || char_line=="Plus" || char_line=="+"){
	char_flag = 1;
    } else {
	throw new Error(text("MarchAlongC: invalid input for char. line"));
    } // end if
    //
    int dir_flag;
    if(dir=="down" || dir=="Down"){
	dir_flag = 1;
    } else if(dir=="up" || dir=="Up"){
	dir_flag = -1;
    } else {
	throw new Error(text("MarchAlongC: invalid input for direction to march"));
    } // end if
    //
    int node1,node2,node4;
    if(char_flag == -1){
	node1 = new_first;
	node2 = old_first;
	nodemarchlist ~= node1;
    } else if(char_flag == 1){
	node1 = old_first;
	node2 = new_first;
	nodemarchlist ~= node2;
    } // end if
    //
    int more = YES;
    do{
	node4 = InteriorNode(node1,node2,-1); //creates new solution node
	nodemarchlist ~= node4;
	if(char_flag == -1){
	    node1 = node4;
	    if(dir_flag == 1){
		node2 = GetNodeData(node2).CMD;
	    } else if(dir_flag == -1){
		node2 = GetNodeData(node2).CMU;
	    } // end if
	} else if(char_flag == 1){
	    node2 = node4;
	    if(dir_flag == 1){
		node1 = GetNodeData(node1).CPD;
	    } else if(dir_flag == -1){node1 = GetNodeData(node1).CPU;
	    } // end if
	} // end if
	if(node2 == NO_NODE || node1 == NO_NODE){more = NO;}
    } while(more == YES);
    //
    return node4;
} // end MarchAlongC()

/**
 * Purpose: Insert a node (node4) in between two initial nodes
 * (node1 and node2).
 * If node1 and node2 are adjacent nodes along a characteristic
 * line, node4 will be connected in between.
 * Input:
 *   node1: index of initial point 1
 *   node2: index of initial point 2
 *   node4: index of solution point (may have a value of -1)
 * If -1 is specified as the index for node4, a new node will be
 * created for the solution point
 *   alpha: fraction that node4 is like node2;
 *   n4.value = alpha * n2.value + (1-alpha) * n1.value
 *   fs: flow state to work in
 * Output:
 * Returns the index of the solution point
 */
int InsertNode(int node1,int node2,int node4,double alpha=0.5,int fs=-1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    //first, assign the Node data
    if(node1 == node2){
	throw new Error(text("InsertNode: node1 & node2 have same index"));
    } // end if
    NodeData n1 = GetNodeData(node1);
    if (n1 is null){
	throw new Error(text("InsertNode: node1 doesn't exist, id: ",node1));
    } // end if
    NodeData n2 = GetNodeData(node2);
    if (n2 is null){
	throw new Error(text("InsertNode: node2 doesn't exist, id: ",node2));
    } // end if
    if (node4==-1){
	node4 = CreateNode(node4);
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("InsertNode: couldn't create node4, id: ",node4));
	} else {
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    // Enforce a 0.0..1.0 range for alpha
    if(alpha > 1.0){alpha = 1.0;}
    if(alpha < 0.0){alpha = 0.0;}
    // Linearly Interpolate all node properties.
    SetNodeData(node4,"X",(1.0 - alpha) * n1.x + alpha * n2.x);
    SetNodeData(node4,"Y",(1.0 - alpha) * n1.y + alpha * n2.y);
    SetNodeData(node4,"P",(1.0 - alpha) * n1.P + alpha * n2.P);
    SetNodeData(node4,"rho",(1.0 - alpha) * n1.rho + alpha * n2.rho);
    SetNodeData(node4,"T",(1.0 - alpha) * n1.T + alpha * n2.T);
    SetNodeData(node4,"V",(1.0 - alpha) * n1.V + alpha * n2.V);
    SetNodeData(node4,"theta",(1.0 - alpha) * n1.theta + alpha * n2.theta);
    // Connect into the mesh only if nodes 1 and 2 are adjacent.
    if(n1.CPD == n2.CPU){
	SetNodeData(node4,"CPlusUp",node1);
	SetNodeData(node1,"CPlusDown",node4);
	SetNodeData(node2,"CPlusUp",node4);
	SetNodeData(node4,"CPlusDown",node2);
    } else if(n1.CPU == n2.CPD){
	SetNodeData(node4,"CPlusUp",node2);
	SetNodeData(node2,"CPlusDown",node4);
	SetNodeData(node1,"CPlusUp",node4);
	SetNodeData(node4,"CPlusDown",node1);
    } else if(n1.CMD == n2.CMU){
	SetNodeData(node4,"CMinusUp",node1);
	SetNodeData(node1,"CMinusDown",node4);
	SetNodeData(node2,"CMinusUp",node4);
	SetNodeData(node4,"CMinusDown",node2);
    } else if(n1.CMU == n2.CMD){
	SetNodeData(node4,"CMinusUp",node2);
	SetNodeData(node2,"CMinusDown",node4);
	SetNodeData(node1,"CMinusUp",node4);
	SetNodeData(node4,"CMinusDown",node1);
    } else if(n1.C0D == n2.C0U){
	SetNodeData(node4,"CZeroUp",node1);
	SetNodeData(node1,"CZeroDown",node4);
	SetNodeData(node2,"CZeroUp",node4);
	SetNodeData(node4,"CZeroDown",node2);
    } else if(n1.C0U == n2.C0D){
	SetNodeData(node4,"CZeroUp",node2);
	SetNodeData(node2,"CZeroDown",node4);
	SetNodeData(node1,"CZeroUp",node4);
	SetNodeData(node4,"CZeroDown",node1);
    } // end if
    // Assuming a successful calculation, return the index of the solution node.
    return node4;
} // end InsertNode()

/**
 * Purpose: Calculate a new streamline node, extending the streamline
 * by length dL 
 * Input:  
 *   node0: index of initial point on the streamline 
 *   node4: index of solution point (may have a value of -1) 
 *   dL: step-size along streamline;
 *     A positive value will step downstream while a negative
 *     value will step upstream. 
 * If -1 is specified as the index for node4, a new node will be
 * created for the solution point. 
 *   fs: flow state to work in
 * Output :  
 * Returns the index of the solution point or a value of -1
 * if there has been a failure.  One possible failure is that
 * there are no nodes close enough to include in the interpolation phase. 
 */
int StepStreamNode(int node0,int node4,double dL,int fs=-1)
{	
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    //
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //first, assign the Node data
    NodeData n0 = GetNodeData(node0);
    if (n0 is null){
	throw new Error(text("StepStreamNode: node0 doesn't exist, id: ",node0));
    } // end if
    int node4_created=NO;
    if (node4 == -1){
	node4 = CreateNode(-1);
	node4_created = YES;
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("StepStreamNode: couldn't create node4, id: ",node4));
	} else {
	    node4_created = YES;
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    //initial guess at solution point, will be very wrong but used to check convergence
    n4.x = n0.x + dL * cos(n0.theta); 
    n4.y = n0.y + dL * sin(n0.theta); 
    n4.theta = n0.theta;
    n4.V = n0.V;
    double Rad = 0.9 * dL; // Radius of influence.
    double mu = 2.0; // smoothing parameter for Shepard interpolation.
    int[] near_node_array;
    NodeData[] near_node_data;
    int near_node_count = FindNodesNear(n4.pos,Rad,near_node_array,MAX_NEAR);
    foreach(int j; near_node_array){
	if(j != node4 && j != node0){
	    NodeData np = GetNodeData(j);
	    near_node_data ~= np;
	} // end if
    } // end foreach
    near_node_count = to!int(near_node_data.length);
    if(near_node_count == 0){
	if(node4_created==YES){
	    DeleteNode(node4);
	} // end if
	return MOC_ERROR;
    } // end if
    // Compute the solution point position and flow properties.
    int iteration_count=0;
    double change_in_position;
    Vector3 old_pos;
    do {
	++iteration_count;
	old_pos = n4.pos;
	// Interpolate the flow properties by Shepard interpolation over the near-by nodes.
	double[] Xi, w;
	double sum_Xi = 0.0;
	n4.V = 0.0;
	n4.theta = 0.0;
	foreach(j; 0 .. near_node_count){
	    NodeData n = near_node_data[j];
	    double r = abs(n4.pos - n.pos);
	    Xi ~= pow((1.0 - r/Rad),mu);
	    sum_Xi += Xi[j];
	} // end foreach
	foreach(j; 0 .. near_node_count){
	    w ~= Xi[j] / sum_Xi;
	    NodeData n = near_node_data[j];
	    n4.V += w[j] * n.V ;
	    n4.theta += w[j] * n.theta ;
	} // end foreach
	// Locate solution point by using average slope.
	double sinCzero = 0.5 * ( sin(n0.theta) + sin(n4.theta) );
	double cosCzero = 0.5 * ( cos(n0.theta) + cos(n4.theta) );
	n4.x = n0.x + cosCzero * dL;
	n4.y = n0.y + sinCzero * dL;
	change_in_position = abs(n4.pos - old_pos);
    }while(iteration_count<MAX_ITERATION && change_in_position>POS_TOL);
    // Save the solution-point properties and connect the node into the characteristic mesh.
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"theta",n4.theta);
    SetNodeData(node4,"V",n4.V);
    SetNodeData(node4,"T",n0.T0 - n4.V^^2*(g-1.)/2./g/R);
    SetNodeData(node4,"P",n0.P0 / p0_p(n4.Mach,g));
    SetNodeData(node4,"rho",n4.P/R/n4.T);
    if (n4.x>n0.x){
	SetNodeData(node4,"CZeroUp",node0);
	SetNodeData(node0,"CZeroDown",node4);
    } else {
	SetNodeData(node4,"CZeroDown",node0);
	SetNodeData(node0,"CZeroUp",node4);
    } // end if
    // Assuming a successful calculation, return the index of the solution node.
    return node4;
} // end StepStreamNode()

/**
 * Purpose: Locate a new node at coordinates (x,y), interpolating
 * the node's properties from other near-by nodes. 
 * Input: 
 *   point: coordinates of the new node 
 *   Rad : radius-of-influence for the Shepard interpolation 
 *   node4 : index of solution point (may have a value of -1) 
 * If -1 is specified as the index for node4, a new node will be
 * created for the solution point. 
 * Output :  
 * Returns the index of the solution point or a value of -1
 * if there has been a failure.  One possible failure is that
 * there are no nodes close enough to include in the interpolation. 
 */
int InterpolateNode(Vector3 point,double Rad,int node4,int fs=-1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    //
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //first, assign the Node data
    int node4_created=NO;
    if (node4 == -1){
	node4 = CreateNode(-1);
	node4_created = YES;
    } // end if
    NodeData n4 = GetNodeData(node4);
    if (n4 is null){
	node4 = CreateNode(node4);
	if(node4 == -1){
	    throw new Error(text("InterpolateNode: couldn't create node4, id: ",node4));
	} else {
	    node4_created = YES;
	    n4 = GetNodeData(node4);
	} // end if
    } // end if
    double mu = 2.0; // smoothing parameter for Shepard interpolation.
    int[] near_node_array;
    NodeData[] near_node_data;
    int near_node_count = FindNodesNear(point,Rad,near_node_array,MAX_NEAR);
    foreach(int j; near_node_array){
	if(j != node4){
	    NodeData np = GetNodeData(j);
	    near_node_data ~= np;
	} // end if
    } // end foreach
    near_node_count = to!int(near_node_data.length);
    if(near_node_count == 0){
	if(node4_created==YES){
	    DeleteNode(node4);
	} // end if
	return MOC_ERROR;
    } // end if
    // The location of the solution point is given.
    n4.pos = point;
    //Interpolate the flow properties by Shepard interpolation over the near-by nodes.
    double[] Xi, w;
    double sum_Xi = 0.0;
    n4.V = 0.0;
    n4.theta = 0.0;n4.T = 0.0;n4.P = 0.0;n4.rho = 0.0;
    foreach(j; 0 .. near_node_count){
	NodeData n = near_node_data[j];
	double r = abs(n4.pos - n.pos);
	Xi ~= pow((1.0 - r/Rad),mu);
	sum_Xi += Xi[j];
    } // end foreach
    foreach(j; 0 .. near_node_count){
	w ~= Xi[j] / sum_Xi;
	NodeData n = near_node_data[j];
	n4.V += w[j] * n.V ;
	n4.theta += w[j] * n.theta ;
    } // end foreach
    // Save the solution-point properties and connect the node into the characteristic mesh.
    SetNodeData(node4,"X",n4.x);
    SetNodeData(node4,"Y",n4.y);
    SetNodeData(node4,"theta",n4.theta);
    SetNodeData(node4,"V",n4.V);
    SetNodeData(node4,"T",n4.T);
    SetNodeData(node4,"P",n4.P);
    SetNodeData(node4,"rho",n4.rho);
    // Assuming a successful calculation, return the index of the solution node.
    return node4;
} // end InterpolateNode

/**
 * Purpose: Creates an oblique shock at a node by creating a new
 * flowstate at this node using the oblique shock relations
 * Requires either the shock angle (beta) or flow turning angle (theta)
 * to be known
 * Input: 
 *   node: index of node 
 *   angle : string of known angle == "beta" || "theta"
 *   value : value of the known angle
 *   fs: flow state for the post-shock data (default = 1)
 */
void MakeOblShkNode(int node,string angle,double value,int fs=1)
{
    if(fs == 0){
	fs = 1;
    } // end if
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    int axiflag=GetAxiFlag();
    //
    NodeData n = GetNodeData(node);
    if (n is null){
	throw new Error(text("MakeOblShkNode: node doesn't exist, id: ",node));
    } // end if
    int addstates = fs-GetNumberFlowStates(node)+1;
    if(addstates > 0){
	AddFlowState(node,addstates);
    } // end if
    double beta,theta;
    if(angle=="beta" || angle=="Beta"){
	beta = value;
	theta = theta_obl(n.F[fs-1].Mach,beta,g);
    } else if(angle=="theta" || angle=="Theta"){
	theta = value;
	beta = beta_obl(n.F[fs-1].Mach,theta,g);
    } else {
	throw new Error(text("MakeOblShkNode: select beta or theta for angle input: ",angle));
    }// end if
    if(axiflag==0){
	SetNodeData(node,"T",n.F[fs-1].T * T2_T1_obl(n.F[fs-1].Mach,beta,g));
	SetNodeData(node,"P",n.F[fs-1].P * p2_p1_obl(n.F[fs-1].Mach,beta,g));
	SetNodeData(node,"rho",n.F[fs-1].rho * r2_r1_obl(n.F[fs-1].Mach,beta,g));
	SetNodeData(node,"theta",n.F[fs-1].theta + theta);
	SetNodeData(node,"V",n.F[fs-1].V * V2_V1_obl(n.F[fs-1].Mach,beta,theta,g));	
    } else if (axiflag==1){
	//TODO add in taylor mcnichol solutions
	SetNodeData(node,"T",n.F[fs-1].T * T2_T1_obl(n.F[fs-1].Mach,beta,g));
	SetNodeData(node,"P",n.F[fs-1].P * p2_p1_obl(n.F[fs-1].Mach,beta,g));
	SetNodeData(node,"rho",n.F[fs-1].rho * r2_r1_obl(n.F[fs-1].Mach,beta,g));
	SetNodeData(node,"theta",n.F[fs-1].theta + theta);
	SetNodeData(node,"V",n.F[fs-1].V * V2_V1_obl(n.F[fs-1].Mach,beta,theta,g));
    } // end if
} // end MakeOblShkNode()

// end of src/moc/unitproc.d

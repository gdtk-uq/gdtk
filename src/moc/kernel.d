/**
 * file: kernel.d
 * location: ./src/moc
 * author: Momar Hughes
 * description: This set of functions forms the computational 
 * and data storage kernel for the IMOC_D program.
 * version: 	
 *		14 Mar 2015: initial port
 *		10 May 2015: use existing classes in CFCFD3 - Vector3
 *		16 May 2015: FindNodesNear now finds closest nodes instead 
 * 					 of just moving through NodeList until done
 *		13 Sep 2015: new node classes NodeData and FlowState
 */
module kernel;

import std.stdio, std.math, std.conv, std.string, std.file;
import gasdynamic;
import gas.gas_model;
import geom;

// Definition of global constants
enum int YES = 1;
enum int NO = 0;
enum int MOC_OK = 0;
enum int MOC_ERROR = -1;
enum int NO_NODE = -1;
//
enum int MAX_NODES = 1000;
enum int MAX_NEAR = 30; //max no. nodes to search for nearby
// definition of default gas properties
double g = 1.4; // ratio of specific heats
double R = 287.0; // gas constant J/mol/K
int AxiFlag = NO; // 1=axisymmetric geometry,0=2D planar
// storage array for nodes
NodeData[MAX_NODES] Node; //array that stores all node data
int NumberOfNodes = 0; //Count no. of valid nodes
int FS = 0; // current flow state to work in

/**
 * FlowState stores gas-state data along with a velocity
 * and flow angle for points in the flow
 */
class FlowState
{
    GasState gas; // gas state
    // [TODO] Momar, I prefer the default D initialization of nan for doubles. PJ
    // [TODO] If they get used uninitialized, it' clear because the nans propagate.
    double theta=0.0; // flow angle, radians
    double V=0.0; // velocity, m/s
    //
    @property ref double rho(){return gas.rho;}
    @property ref double P(){return gas.p;}
    @property ref double T(){return gas.T[0];}
    @property double a()
    {
    	return sqrt(g*R*T);
    } // end a()
    @property double Mach()
    {
    	return V/a;
    } // end Mach()
    @property double Nu()
    {
    	return PM1(Mach,g);
    } // end Nu()
    @property double T0()
    {
    	return T * T0_T(Mach,g);
    } // end T0()
    @property double P0()
    {
    	return P * p0_p(Mach,g);
    } // end P0()
    //	
    this()
    {
	this.gas = new GasState(1, 1); // PJ single-species gas, 1 internal energy mode
	// [TODO] We should construct this gas state with the knowledge of the gas model.
	this.gas.p = 0.0;
	this.gas.rho = 0.0;
	// this.gas.T ~= [0.0]; // [TODO] don't like this; it makes the array bigger.
    } // end this()
    //
    this(GasState gas,double theta,double V)
    {
	this.gas = new GasState(gas);
	this.theta = theta;
	this.V = V;
    } // end this()
    //
    this(Vector3 V)
    {
	this.V = sqrt(V.x^^2+V.y^^2);
	this.theta = atan(V.y/V.x);
    } // end this()
    //
    this(FlowState other)
    {
	this.gas = new GasState(other.gas);
	this.theta = other.theta;
	this.V = other.V;
    } // end this()
    //
    FlowState opAssign(double theta,double V)
    {
	this.theta = theta;
	this.V = V;
	return this;
    } // end opAssign()
    //
    FlowState opAssign(GasState gas)
    {
	this.gas = new GasState(gas);
	return this;
    } // end opAssign()
    //
    FlowState dup()
    {
    	return new FlowState(gas,theta,V);
    } // end dup()
    //
} // end class FlowState

/**
 * NodeData stores the data for nodes which have one position
 * but may have multiple flow-states
 */
class NodeData 
{
public:
    Vector3 pos = Vector3(0.0,0.0); // node coordinates
    FlowState[] F; // array of flowstates at node
    int CPU=NO_NODE;     // upstream node along C+ characteristic
    int CPD=NO_NODE;   // downstream node along C+ characteristic
    int CMU=NO_NODE;    // upstream node along C- characteristic
    int CMD=NO_NODE;  // downstream node along C- characteristic
    int C0U=NO_NODE;     // upstream along streamline
    int C0D=NO_NODE;   // downstream along streamline
    // quick access to position vector
    @property ref double x(){return pos._p[0];}
    @property ref double y(){return pos._p[1];}
    // quick access to first flowstate
    @property ref FlowState F0(){return F[0];}    
    @property ref FlowState F1(){return F[1];}
    @property ref FlowState F2(){return F[2];}
    //
    @property ref double rho(){return F[FS].gas.rho;}
    @property ref double P(){return F[FS].gas.p;}
    @property ref double T(){return F[FS].gas.T[0];}
    @property ref double V(){return F[FS].V;}
    @property ref double theta(){return F[FS].theta;}
    @property double a(){return F[FS].a;}
    @property double Mach(){return F[FS].Mach;}
    @property double Nu(){return F[FS].Nu;}
    @property double T0(){return F[FS].T0;}
    @property double P0(){return F[FS].P0;}
    //
    void AddFlowState()
    {
    	F ~= new FlowState;
    } // end AddFlowState()
    void DeleteFlowState()
    {
    	F.length = F.length-1;
    } // end DeleteFlowState()
    //
    this()
    {
    	this.F ~= new FlowState;
    } // end this
    //
    this(Vector3 pos,FlowState[] F,int CPU,int CPD,int CMU,int CMD,int C0U,int C0D)
    {
    	this.pos = pos.dup;
	this.F = F.dup;
	this.CPU = CPU;
	this.CPD = CPD;
	this.CMU = CMU;
	this.CMD = CMD;
	this.C0U = C0U;
	this.C0D = C0D;
    } // end this()
    //
    this(Vector3 pos,FlowState F,int CPU,int CPD,int CMU,int CMD,int C0U,int C0D)
    {
	this.pos = pos.dup;
	this.F[FS] = F.dup;
	this.CPU = CPU;
	this.CPD = CPD;
	this.CMU = CMU;
	this.CMD = CMD;
	this.C0U = C0U;
	this.C0D = C0D;
    } // end this()
    //
    this(Vector3 pos,FlowState[] F)
    {
	this.pos = pos.dup;
	this.F = F.dup;
    } // end this()
    //
    this(Vector3 pos,FlowState F,int fs=0)
    {
	this.pos = pos;
	this.F[fs] = F.dup;
    } // end this()
    //
    NodeData opAssign(Vector3 pos)
    {
	this.pos = pos;
	return this;
    } // end opAssign()
    //
    NodeData opAssign(FlowState F,int state=0)
    {
	this.F[state] = F;
	return this;
    } // end opAssign()
    //
    NodeData dup()
    {
    	return new NodeData(pos,F,CPU,CPD,CMU,CMD,C0U,C0D);
    }
    //
    @property double beta(int state=FS)
    {
	if(state<1){
	    throw new Error(text("node does not have enough flowstates to calculate beta"));
	} else {
	    double M = F[state-1].Mach;
	    double theta = F[state].theta - F[state-1].theta;
	    return beta_obl(M,theta,g);
	} // end if
    } // end beta()
    //
}//end class NodeData

/**
 * Sets the ratio of specific heats (assumed constant 
 * throughout flow field)
 * Input:
 *    value: ratio of specific heats
 */
void SetRatioOfSpecificHeats(double value)
{
    g = value;
    if(g < 0.0){
	throw new Error(text("SetGasConstant: invalid value : ",value));
    } // end if
} //end SetRatioOfSpecificHeats()

/**
 * Returns the current ratio of specific heats
 */
double GetRatioOfSpecificHeats()
{
    return g;
} //end GetRatioOfSpecificHeats()

/**
 * Sets the gas constant (assumed constant throughout
 * flow field)
 * Input:
 *    value: gas constant (J/mol/Kg)
 */
void SetGasConstant(double value)
{
    R = value;
    if(R < 0.0){
	throw new Error(text("SetGasConstant: invalid value : ",value));
    } // end if
} //end SetGasConstant()

/**
 * Returns the current gas constant (J/mol/Kg)
 */
double GetGasConstant()
{
    return R;
}//end GetGasConstant()

/**
 * Sets the axisymmetric analysis flag
 * Input:
 *    value: 1: axisymmetric, 0: 2D planar
 */
void SetAxiFlag(int value)
{
    AxiFlag = value;
    if(AxiFlag != 0 && AxiFlag != 1){
	throw new Error(text("SetAxiFlag: invalid value : ",value));
    } // end if
} //end SetAxiFlag

/**
 * Returns the current axisymmetric analysis flag
 * 1: axisymmetric, 0: 2D planar
 */
int GetAxiFlag()
{
    return AxiFlag;
} //end GetAxiFlag

/**
 * Returns the total current number of valid nodes
 */
int GetNumberOfNodes()
{
    return NumberOfNodes;
} //end GetNumberOfNodes

void SetFlowState(int value)
{
    if(value < 0){
	throw new Error(text("SetFlowState: invalid value : ",value));
    } // end if
    FS = value;
} // end SetFlowState()

int GetFlowState()
{
    return FS;
} // end GetFlowState()

/**
 * Checks that a node is valid; is not null
 * Input:
 *    id: node id
 */
int ValidNode(int id)
{
    if(id>=0 && id <= MAX_NODES && Node[id] !is null){
        return YES;
    } else {
        return NO;
    } // end if
} //end ValidNode

/**
 * Provides access to the NodeData object at a requested id
 * Input:
 *    id: node id
 */
NodeData GetNodeData(int id)
{
    if(ValidNode(id)==YES){
        return Node[id];
    } else {
        //throw new Error(text("GetNodeData: invalid node, id: ",id));
        return null;
    } // end if
} // end GetNodeData()

/**
 * Prints important node properties to screen
 * Input:
 *   id: node id
 */
void WriteNodeData(int id)
{
    if(ValidNode(id)==YES){
	auto n = GetNodeData(id);
	writefln("Node %d: x=%f,y=%f,C0D=%d,C0U=%d,CMD=%d,CMU=%d,CPD=%d,CPU=%d" ,id,n.x,n.y,n.C0D,n.C0U,n.CMD,n.CMU,n.CPD,n.CPU);
	foreach(j;0 .. n.F.length){
	    if(j>0){writefln("beta=%f",n.beta*180./PI);}
	    writefln("  state %d: P0=%f,T0=%f,p=%f,T=%f,rho=%f,theta=%f,V=%f,Mach=%f,Nu=%f" ,j,n.F[j].P0,n.F[j].T0,n.F[j].P,n.F[j].T,n.F[j].rho,n.F[j].theta*180.0/PI,n.F[j].V,n.F[j].Mach,n.F[j].Nu*180.0/PI);
	}
    } else {
	writefln("Node %d is not a valid node",id);
    } // end if
}//end WriteNodeData

/**
 * Creates a new node at the selected if
 * Input:
 *   id < 0: Select the next available location
 *   id >= 0: Use a specific location
 * Output:
 * Returns: index of the created node
 *          -1 on failure. 
 */
int CreateNode(int id=-1,int fs=-1)
{
    if(fs != -1){
	SetFlowState(fs);
    } // end if
    fs = GetFlowState();
    int foundSpace;
    if(id>=MAX_NODES){
	throw new Error(text("CreateNode: invalid node id: ",id));
    } // end if
    if(id<0){ //search for next vacant space in Node
	foundSpace=NO;
	foreach(i;0 .. MAX_NODES){
	    if(Node[i] is null){
		foundSpace = YES;
		id = i; break;
	    } // end if
	} // end foreach
	if(foundSpace == NO){ //Node array is already full
            throw new Error(text("ERROR - CreateNode: Maximum number of nodes reached"));
	} // end if
    } else { //Node id has been explicitly specified
	if(Node[id] !is null){ //clear id address if currently full
	    //writefln("CreateNode: Warning; Node %d data deleted and replaced",id);
	    DeleteNode(id);
	} // end if
    } // end if
    Node[id] = new NodeData;
    AddFlowState(id,fs);
    if(Node[id] !is null){
	++NumberOfNodes;
	return id;
    } else { //CreateNode failed, address could not be allocated
	throw new Error(text("CreateNode: Memory allocation failed"));
    } // end if
} // end CreateNode

/**
 * Destroys the memory allocated to a node and relinks 
 * any connected nodes
 * Input:
 *    id: node id
 */ 
void DeleteNode(int id)
{
    if(ValidNode(id)==YES){
	int idCPU=Node[id].CPU, idCPD=Node[id].CPD;
	if(idCPU != NO_NODE && Node[idCPU] !is null){
	    Node[idCPU].CPD=idCPD;
	} // end if
	if(idCPD != NO_NODE && Node[idCPD] !is null){
	    Node[idCPD].CPU=idCPU;
	} // end if
	int idCMU=Node[id].CMU, idCMD=Node[id].CMD;
	if(idCMU != NO_NODE && Node[idCMU] !is null){
	    Node[idCMU].CMD=idCMD;
	} // end if
	if(idCMD != NO_NODE && Node[idCMD] !is null){
	    Node[idCMD].CPU=idCMU;
	} // end if
	int idC0U=Node[id].C0U, idC0D=Node[id].C0D;
	if(idC0U != NO_NODE && Node[idC0U] !is null){
	    Node[idC0U].C0D=idC0D;
	} // end if
	if(idC0D != NO_NODE && Node[idC0D] !is null){
	    Node[idC0D].C0U=idC0U;
	} // end if
	Node[id] = null;
	--NumberOfNodes;  
    } // end if
}//end DeleteNode()

/**
 * Returns the number of flowstates present at a node
 * Input:
 *    id: node id
 */
int GetNumberFlowStates(int id)
{
    if(ValidNode(id)==YES){
        return to!int(Node[id].F.length);
    } else {
        return 0;
    } // end if
} // end GetNumberFlowStates()

/**
 * Creates a new flowstate at a node
 * Input:
 *    id: node id
 */
void AddFlowState(int id,int states=1)
{
    if(ValidNode(id)==YES){
    	foreach(i;0 .. states){
	    Node[id].AddFlowState;
        } // end foreach
    } // end if
} // end AddFlowState()

/**
 * Destroys the last flowstate at a node
 * Input:
 *    id: node id
 */
void DeleteFlowState(int id)
{
    if(ValidNode(id)==YES){
        Node[id].DeleteFlowState;
    }
} // end DeleteFlowState()

/**
 * Set the value for a particular node variable
 * Inputs:
 *   variable: string specifying which variable to set
 *   dvalue: pointer to a string specifying value for the variable
 *   id: index of the node
 *   state: index of flowstate set value for
 */ 
void SetNodeData(int id,string variable,double dvalue,int fs=-1)
{
    if(fs == -1){
	fs = GetFlowState();
    } // end if
    if(ValidNode(id)==YES){
	//find, then change, the correct Node variable
        int ivalue = to!int(dvalue); //int type of value
        while(Node[id].F.length-1 < fs){Node[id].AddFlowState;}
        if(variable=="X"){
	    Node[id].x = dvalue;
	} else if(variable=="Y"){
	    Node[id].y = dvalue;
	} else if(variable=="P"||variable=="p"){
	    Node[id].F[fs].P = dvalue;
	} else if(variable=="T"){
	    Node[id].F[fs].T = dvalue;
	} else if(variable=="rho"){
	    Node[id].F[fs].rho = dvalue;
	} else if(variable=="theta"||variable=="Theta"){
	    Node[id].F[fs].theta = dvalue;
	} else if(variable=="V"||variable=="Vel"){
	    Node[id].F[fs].V = dvalue;
	} else if(variable=="CPlusUp"||variable=="CPU"){
	    if(ivalue==id){
		throw new Error(text("SetNodeData: Node cannot reference itself, id: "));
	    } // end if
	    Node[id].CPU = ivalue;
	} else if(variable=="CPlusDown"||variable=="CPD"){
	    if(ivalue==id){
		throw new Error(text("SetNodeData: Node cannot reference itself, id: "));
	    } // end if
	    Node[id].CPD = ivalue;
	} else if(variable=="CMinusUp"||variable=="CMU"){
	    if(ivalue==id){
		throw new Error(text("SetNodeData: Node cannot reference itself, id: "));
	    } // end if
	    Node[id].CMU = ivalue;
	} else if(variable=="CMinusDown"||variable=="CMD"){
	    if(ivalue==id){
		throw new Error(text("SetNodeData: Node cannot reference itself, id: "));
	    } // end if
	    Node[id].CMD = ivalue;
	} else if(variable=="CZeroUp"||variable=="C0U"){
	    if(ivalue==id){
		throw new Error(text("SetNodeData: Node cannot reference itself, id: "));
	    } // end if
	    Node[id].C0U = ivalue;
	} else if(variable=="CZeroDown"||variable=="C0D"){
	    if(ivalue==id){
		throw new Error(text("SetNodeData: Node cannot reference itself, id: "));
	    } // end if
	    Node[id].C0D = ivalue;
	} else {
	    throw new Error(text("SetNodeData: invalid variable string: ",variable));
	} // end if
    } else {
        throw new Error(text("SetNodeData: invalid node, id: ",id));
    } // end if
} // end SetNodeData()

/**
 * Search for the next Node that exists after the specified starting position
 * Input: 
 *   idStart : Node index to start searching from
 * 				-1 searches from and including the first node index
 * Output :
 *   id : next Node id
 */
int GetNextNodeId(int idStart)
{
    if(idStart<-1){
        idStart = -1; //correct invalid id
    } // end if
    if(idStart >=(MAX_NODES -1)){
        writeln("GetNextnodeId: started search from last node");
        return MOC_ERROR;
    } // end if
    foreach(id;idStart+1 .. MAX_NODES){
	if(Node[id] !is null){
	    return id;
	} // end if
    } // end foreach
    return MOC_ERROR; // no next Node found
} // end GetNextNodeId()

/**
 * Search for the nearest nodes to a particular point
 * Input: 
 *   point: coordinates to seach around
 *   tol: distance tolerance within which to search
 *		  a negative tolerance will return the single closest node
 *   idNearArray: array to save indices of nearby nodes
 *   maxCount: maximum number of nearby nodes to save
 * Output :
 *   Returns number of nearby nodes found
 */
int FindNodesNear(Vector3 point, double tol, out int[] idNearArray, int maxCount=MAX_NEAR)
{
    double distNear=10.0e6; //something v large;
    int idNear=-1, idcount=0, nodeCount=0;
    double[] dist_list;
    int[] id_list;
    foreach(id;0 .. MAX_NODES){
	if(Node[id] !is null){
	    NodeData n = GetNodeData(id);
	    double dist = abs(n.pos - point);
	    dist_list ~= dist;
	    id_list ~= id;
	    ++idcount;
	    if(dist<distNear){ //update nearest Node found
		distNear = dist;
		idNear = id;
	    } // end if
	} // end if
    } // end foreach
    //
    if(tol<=0.0){//find nearest Node
	idNearArray ~= idNear;
	nodeCount = 1;
	return nodeCount;
    } else {
	foreach(j;0 .. maxCount){
	    distNear = 10e6;
	    foreach(id;0 .. idcount){
		if(dist_list[id]<=distNear){ //update nearest Node found
		    distNear = dist_list[id];
		    idNear = id;
		} // end if
	    } // end foreach
	    if(distNear<=tol){
		idNearArray ~= idNear;
		++nodeCount;
		dist_list[idNear]=10.0e6;
	    } // end if
	} // end foreach
    } // end if	
    return nodeCount;
} // end FindNodesNear()		

/**
 * Creates a string summarising the nearest nodes to a particular point
 * Input: 
 *   point: coordinates to seach around
 *   tol: distance tolerance within which to search
 *		  a negative tolerance will return the single closest node
 *   maxCount: maximum number of nearby nodes to save
 * Output :
 *   Returns a formatted string of the nearby node indices
 */
string ListNodesNear(Vector3 point,double tol,int maxCount=MAX_NEAR)
{
    int[] id_array;
    int nodeCount = FindNodesNear(point,tol,id_array,maxCount);
    //
    string idString;
    foreach(id;0 .. nodeCount){
	idString ~= to!string(id_array[id]); //CHECK
	if(id == nodeCount-1){break;}
	idString ~= ", ";
    } // end foreach
    return idString;
}//end ListNodesNear()

/**
 * Checks whether two nodes are on the same char. line or streamline
 * Input: 
 *   node1: index of first node
 *   node2: index of the second node
 *   line_flag: -1/1/0 for C- char./C+ char./streamline
 * Output :
 *   Returns 1/0 if nodes are/are not on the same line
 */
int CheckSameLine(int node1,int node2,int line_flag)
{
    if(Node[node1] is null){
	throw new Error(text("CheckSameLine: invalid node1, id: ",node1));
    } // end if
    if(Node[node2] is null){
	throw new Error(text("CheckSameLine: invalid node2, id: ",node2));
    } // end if
    int upstream,downstream,node1_new=node1;
    
    while(upstream != -1 && upstream != node2){
	if(line_flag==0){upstream = Node[node1_new].C0U;}
	if(line_flag==-1){upstream = Node[node1_new].CMU;}
	if(line_flag==1){upstream = Node[node1_new].CPU;}
	if(upstream==-1){break;}
	node1_new = upstream;
	writeln("node1_new",node1_new,"upstream",upstream);
    } 
    if(node1_new==node2){
	return YES;
    } // end if
    node1_new=node1;
    while(downstream != -1 && downstream != node2){
	if(line_flag==0){downstream = Node[node1_new].C0D;}
	if(line_flag==-1){downstream = Node[node1_new].CMD;}
	if(line_flag==1){downstream = Node[node1_new].CPD;}
	if(downstream==-1){break;}
	node1_new = downstream;
    }
    if(node1_new==node2){
	return YES;
    } else {
	return NO;
    } // end if
} // CheckSameLine()

/**
 * Saves all important node data to a text file that can later be reloaded
 * Input: 
 *   FileName: name of file to save to
 */
void SaveNodes(string FileName)
{
    File fp = File(FileName, "w"); //write only access to FileName
    fp.writeln("# Id state X Y p T rho Theta V CPlusUp CMinusUp CZeroUp CPlusDown CMinusDown CZeroDown");
    //write data for all Nodes as a string"
    foreach(id;0 .. MAX_NODES){
	NodeData n = Node[id];
	if(n is null){continue;}//move on to next Node
        foreach(j;0 .. n.F.length){
            fp.writefln("%d %d %e %e %e %e %e %e %e %d %d %d %d %d %d",
			id,j+1,n.x,n.y,n.F[j].P,n.F[j].T,n.F[j].rho,
			n.F[j].theta,n.F[j].V,n.CPU,n.CMU,n.C0U,n.CPD,n.CMD,n.C0D);
	} // end foreach
    } // end foreach
    fp.close();
}//end SaveNodes()

/**
 * Load all important node data from a text file created using SaveNodes
 * Input: 
 *   FileName: name of file to load from
 */
void LoadNodes(string FileName)
{
    if(exists(FileName) == 0){ //file does not exist
	throw new Error(text("LoadNode: file does not exist: ",FileName));
    } // end if
    if(isFile(FileName) == 0){
	throw new Error(text("LoadNode: not a valid file format: ",FileName));
    } // end if
    //	
    File fp = File(FileName, "r");
    double value;
    string[] variables=["X","Y","p","T","rho","Theta","V","CPU","CMU","C0U","CPD","CMD","C0D"];
    while(!fp.eof()){
	auto line = split(fp.readln());
	if(line.length != 15){continue;}
	int id = to!int(line[0]);
	int fs = to!int(line[1]);
	if(Node[id] is null){CreateNode(id);}
	foreach(int i;0 .. 13){
	    value = to!double(line[i+2]);
	    SetNodeData(id,variables[i],value,fs);
	} // end foreach
    } // end while
} // LoadNodes()

// END KERNEL

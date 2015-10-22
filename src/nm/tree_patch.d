
/**
 * tree_patch.d
 *
 * Contains classes that can be used to build a look-up table (LUT) for a bi-variate function F(x,y). The table can then be searched.
 * The method divides the x,y domain into patches organised in a binomial tree. 
 * Patches are adaptively sized based on the level of error in that patch.
 * The method is based on "Fast Evaluation of Complex Equations of State" by Luke & Collins (2013)
 *
 * Author: Jonathan Ho
 * Date: 17-09-2015
 */
module nm.tree_patch;
import std.stdio;
import std.math;
import std.mathspecial;
import std.algorithm; 
import std.string;
import std.conv;
import std.datetime;


alias F_xy = double function(double, double);
alias F_transform = double[2] function(double, double);

static double[16][16] B_inv = 
	[[1.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
	[-0.8333333333333335, 3.0000000000000000, -1.5000000000000000, 0.3333333333333333, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
	[0.3333333333333334, -1.5000000000000000, 3.0000000000000000, -0.8333333333333333, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
	[0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 1.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
	[-0.8333333333333335, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 3.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -1.5000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.3333333333333333, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
	[0.6944444444444451, -2.5000000000000004, 1.2500000000000002, -0.2777777777777777, -2.5000000000000004, 9.0000000000000000, -4.5000000000000000, 0.9999999999999997, 1.2500000000000002, -4.5000000000000000, 2.2500000000000000, -0.4999999999999997, -0.2777777777777777, 0.9999999999999998, -0.4999999999999997, 0.1111111111111111],
	[-0.2777777777777780, 1.2500000000000002, -2.5000000000000004, 0.6944444444444444, 1.0000000000000002, -4.5000000000000000, 9.0000000000000000, -2.4999999999999996, -0.5000000000000001, 2.2500000000000000, -4.5000000000000000, 1.2499999999999996, 0.1111111111111111, -0.4999999999999998, 0.9999999999999996, -0.2777777777777777],
	[0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.8333333333333335, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 3.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -1.5000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.3333333333333333],
	[0.3333333333333334, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -1.5000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 3.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.8333333333333333, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
	[-0.2777777777777780, 1.0000000000000002, -0.5000000000000001, 0.1111111111111111, 1.2500000000000002, -4.5000000000000000, 2.2500000000000000, -0.4999999999999998, -2.5000000000000004, 9.0000000000000000, -4.5000000000000000, 0.9999999999999996, 0.6944444444444444, -2.4999999999999996, 1.2499999999999996, -0.2777777777777778],
	[0.1111111111111112, -0.5000000000000001, 1.0000000000000002, -0.2777777777777778, -0.5000000000000001, 2.2500000000000000, -4.5000000000000000, 1.2499999999999996, 1.0000000000000002, -4.5000000000000000, 9.0000000000000000, -2.4999999999999991, -0.2777777777777778, 1.2499999999999996, -2.4999999999999991, 0.6944444444444443],
	[0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.3333333333333334, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -1.5000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 3.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.8333333333333333],
	[0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 1.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
	[0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.8333333333333335, 3.0000000000000000, -1.5000000000000000, 0.3333333333333333],
	[0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.3333333333333334, -1.5000000000000000, 3.0000000000000000, -0.8333333333333333],
	[0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 1.0000000000000000]];
class patch{
	double x_lo;
	double x_hi;
	double y_lo;
	double y_hi;
	double[16] bs;//controlpoints
	
	const double x_m(){return (x_lo + x_hi)*0.5;}
	
	const double y_m(){return (y_lo + y_hi)*0.5;}
		
	double aspectRatio(){return (y_hi-y_lo)/(x_hi-x_lo);}
	
	double area(){return (y_hi-y_lo)*(x_hi-x_lo);}
	
	double transformedArea(F_transform F){
		/*
		p1-------p2
		|         |
		|         |
		|         |
		p0-------p3
		*/
		//transforms x, y into some other two-dimensional plane based on F and calculates the area
		//of the quadilateral in this area
		double[2] p0 = F(x_lo,y_lo);
		double[2] p1 = F(x_lo, y_hi);
		double[2] p2 = F(x_hi, y_hi);
		double[2] p3 = F(x_hi, y_lo);
		double[2] A = [p1[0] - p0[0], p1[1] - p0[1]];
		double[2] B = [p2[0] - p0[0], p2[1] - p0[1]];
		double[2] C = [p3[0] - p0[0], p3[1] - p0[1]];
		return 0.5*(fabs(A[0]*B[1]-B[1]*A[0]) + fabs(B[0]*C[1]-B[1]*C[0]));
	}		
	override string toString() const
	{
		return format("This is a patch with x_lo: %s, x_hi %s, y_lo: %s, y_hi: %s", x_lo, x_hi, y_lo, y_hi);
	}
	this(in double x_lo, in double x_hi, in double y_lo, in double y_hi)
	{
		this.x_lo = x_lo;
		this.x_hi = x_hi;
		this.y_lo = y_lo;
		this.y_hi = y_hi;		
	}
	const double interpolateF(double x, double y){
		assert(((x>=x_lo)&&(x<=x_hi)),"x not bounded properly");
		assert(((y>=y_lo)&&(y<=y_hi)),format("y not bounded properly, y: %s is not in the bracket [%s, %s]",y,y_lo,y_hi));
		double u = (x - x_lo)/(x_hi-x_lo);
		double v = (y - y_lo)/(y_hi-y_lo);
		double[4][4][4] bs_k;//now with the k-th dimension
		int i_bs = 0;
		foreach(ref row;bs_k[3]){
			foreach(ref element; row){
				element = bs[i_bs];
				i_bs++;
				}
		}
		for(int k = 2; k != -1; k--){
			for(int i = 0; i != k+1; i++){
				for(int j = 0; j != k+1; j++){
					bs_k[k][j][i] = (1-u)*(1-v)*bs_k[k+1][j][i]+u*(1-v)*bs_k[k+1][j][i+1]+(1-u)*v*bs_k[k+1][j+1][i] + u*v*bs_k[k+1][j+1][i+1];
					assert(!isNaN(bs_k[k][j][i]),format("i: %s, j: %s, k: %s",i,j,k));
				}
			}
		}
		return bs_k[0][0][0];
		}
	
	void getControlPoints(F_xy F)
		{
		double u;
		double v;
		double x;
		double y;
		double[16] fs;
		for(int n = 0; n!=4; n++ ){
			for(int m = 0; m!=4; m++){//parallelisation possible here
				u = m/3.0; v = n/3.0;
				x = (1.0 - u)*this.x_lo + u*this.x_hi;
				y = (1.0 - v)*this.y_lo + v*this.y_hi;
				fs[4*n+m] = F(x,y);
				if (isNaN(fs[4*n+m])) throw new Exception(format("Function calculated NaN at x: %s, y: %s", x,y));
			}
		}
		foreach(i, ref b; this.bs){
			b = 0.0;
			foreach(j, f; fs) b += B_inv[i][j]*f;
			}
		
		}
	double evalError(F_xy F, double npoints = 10)
	{	
		//EVALUATES ERROR ON A MEAN-SQUARED ERROR METHOD
		//ON A GRID OF npoints x npoints
		double error2 = 0;
		double x;
		double y;
		
		for(int i = 0; i != npoints; i++){
			for(int j = 0; j != npoints; j++){
				x = i/npoints*(this.x_hi - this.x_lo) + x_lo;
				y = j/npoints*(this.y_hi - this.y_lo) + y_lo;
				error2 += pow(this.interpolateF(x,y)/F(x,y) - 1,2);
				
				
			}
		}
		return sqrt(error2/npoints/npoints);
	}
	double evalMaxError(F_xy F, double npoints = 20)
	{	
		//EVALUATES MAX ERROR BY SAMPLING npoints x npoints
		double maxError = 0;
		double x;
		double y;
		
		for(int i = 0; i != npoints; i++){
			for(int j = 0; j != npoints; j++){
				x = i/npoints*(this.x_hi - this.x_lo) + x_lo;
				y = j/npoints*(this.y_hi - this.y_lo) + y_lo;
				maxError = max(fabs(this.interpolateF(x,y)/F(x,y) - 1),maxError);
			}
		}
		return maxError;
	}	
	
	patch splitX_L(){
		double x_m = (x_lo + x_hi)/2.0;
		return new patch(x_lo, x_m, y_lo, y_hi);
	}
	patch splitX_R(){
		double x_m = (x_lo + x_hi)/2.0;
		return new patch(x_m, x_hi, y_lo, y_hi);
		}
	patch splitY_L(){
		double y_m = (y_lo + y_hi)/2.0;
		return new patch(x_lo,x_hi,y_lo,y_m);
	}
	patch splitY_R(){
		double y_m = (y_lo + y_hi)/2.0;
		return new patch(x_lo,x_hi,y_m,y_hi);
	}
	
} 


class TreeNode {
	patch nodePatch;
	TreeNode* left;
	TreeNode* right;
	int idx;
	char splitID;
	
	this(patch myPatch){
		this.nodePatch = myPatch;}
	override string toString() const 
	{
		return format("Node with idx: %s, splitID: %s, left_idx: %s, right_idx: %s", idx, splitID, left, right);
	}
	void writeData(string filename = "Tree.dat",char lastFlag = 'N')
	{
		File tFile = File(filename, "a");
		tFile.writeln(idx);
		if (splitID != 'N') {
			tFile.writeln((*left).idx);
			tFile.writeln((*right).idx);
			}
		else{
			tFile.writeln(0);
			tFile.writeln(0);
			}
		tFile.writeln(splitID);
		tFile.writefln("%.16f",nodePatch.x_lo);
		tFile.writefln("%.16f",nodePatch.x_hi);
		tFile.writefln("%.16f",nodePatch.y_lo);
		tFile.writefln("%.16f",nodePatch.y_hi);
		foreach(b; nodePatch.bs) tFile.writef("%.16e ", b);
		if (lastFlag == 'N') tFile.write("\n");
				
	}
}
class Tree {
	TreeNode[] Nodes;//the actual tree
	double x_lo;//overall patch that the Tree covers
	double x_hi;
	double y_lo;
	double y_hi;
	//OPTIONAL FOR STORING THE ORIGINAL BOUNDS THAT THE TREE WAS TRANSFORMED FROM
	double X_min;
	double X_max;
	double Y_min;
	double Y_max;
	double globalAspectRatio;
	double globalMinArea;
	double globalMaxError;
	
	//----------------------------------------------------------
	//constructor which fills the tree with its first node
	this(double x_lo, double x_hi, double y_lo, double y_hi)
	{
		this.globalAspectRatio = (y_hi-y_lo)/(x_hi - x_lo);
		this.globalMinArea = 1e7; 
		this.globalMaxError = 0.001;
		this.x_lo = x_lo;
		this.x_hi = x_hi;
		this.y_lo = y_lo;
		this.y_hi = y_hi;
		Nodes ~= new TreeNode(new patch(x_lo, x_hi, y_lo, y_hi));//append the first patch
		
		}
	
	void grow(F_xy F, F_transform F_t, TreeNode CurrentTreeNode){
	//grows the tree based on the last node
		int idx = to!int(Nodes.length-1);
		//writefln("-------------------------idx: %s -----------------------------", idx);
		CurrentTreeNode.idx = idx;
		if (Nodes.length == 0){throw new Exception("There are no nodes in the tree");}
		CurrentTreeNode.nodePatch.getControlPoints(F);//get control points to evaluate error
		if((CurrentTreeNode.nodePatch.evalError(F) > globalMaxError)&(CurrentTreeNode.nodePatch.transformedArea(F_t) > globalMinArea)){
			if (CurrentTreeNode.nodePatch.aspectRatio() < globalAspectRatio){
				CurrentTreeNode.splitID = 'X';
				this.Nodes ~= new TreeNode(CurrentTreeNode.nodePatch.splitX_L);
				CurrentTreeNode.left = &Nodes[$-1]; 
				this.grow(F,F_t,Nodes[$-1]);
				
				this.Nodes ~= new TreeNode(CurrentTreeNode.nodePatch.splitX_R);
				CurrentTreeNode.right = &Nodes[$-1];
				this.grow(F,F_t,Nodes[$-1]);
			}
			else{
				CurrentTreeNode.splitID = 'Y';
				this.Nodes ~= new TreeNode(CurrentTreeNode.nodePatch.splitY_L);
				CurrentTreeNode.left = &Nodes[$-1];
				this.grow(F,F_t,Nodes[$-1]);
				this.Nodes ~= new TreeNode(CurrentTreeNode.nodePatch.splitY_R);
				CurrentTreeNode.right = &Nodes[$-1];
				this.grow(F,F_t,Nodes[$-1]);
				}
			}
		else{
			CurrentTreeNode.splitID = 'N';
			}
	}
	const const(patch) search(double x, double y){
	//returns the patch to the node in the tree that bounds x and y
	//had to use pointers as it was const
	const(TreeNode)* currentNode = &Nodes[0];//start at first node
		if ((x<currentNode.nodePatch.x_lo)||(x>currentNode.nodePatch.x_hi)){
			throw new Exception(format("x (probably re-parameterized internal energy) not bounded properly by look up table, x: %s is not in the bracket [%s, %s]",x,currentNode.nodePatch.x_lo,currentNode.nodePatch.x_hi));
		}
		if ((y<currentNode.nodePatch.y_lo)||(y>currentNode.nodePatch.y_hi)){
			throw new Exception(format("y (probably re-parametrized density) not bounded properly by look up table, y: %s is not in the bracket [%s, %s]",y,currentNode.nodePatch.y_lo,currentNode.nodePatch.y_hi));
		}
		while((*currentNode).splitID != 'N') {
			if (currentNode.splitID == 'X'){
				if (x < (*currentNode).nodePatch.x_m()) currentNode = (*currentNode).left;
				else currentNode = (*currentNode).right;
				}
			else if (currentNode.splitID == 'Y'){
				if (y < (*currentNode).nodePatch.y_m()) currentNode = (*currentNode).left;
				else currentNode = (*currentNode).right;
				}
			
			} 
		return (*currentNode).nodePatch;
	}
	
	void writeLeaves(){
		//write's co-ordinates of all leaf node patches for plotting in python
		File patchFile = File("patch_xy.dat", "w");
		foreach (node; Nodes){
			if (node.splitID == 'N'){
				patchFile.writeln([node.nodePatch.x_lo, node.nodePatch.x_hi,node.nodePatch.y_lo,node.nodePatch.y_hi]);
				}
			}
		
		}
	void writeLUT(string filename = "Tree.dat"){
		File tFile = File(filename, "w");
		tFile.writeln(globalMaxError);
		tFile.writeln(globalMinArea);
		tFile.writeln(X_min);
		tFile.writeln(X_max);
		tFile.writeln(Y_min);
		tFile.writeln(Y_max);
		tFile.close();
		foreach(i, Node; Nodes) {
			if(i != Nodes.length - 1) Node.writeData(filename);
			else Node.writeData(filename,'L');
		}
		}
}

Tree buildTree_fromFile(string filename = "Tree.dat"){
	File treeFile = File(filename,"r");
	int i = 0;
	string[9] lines;
	int idx;
	int[] lefts;
	int[] rights;
	char splitID;
	double x_lo;
	double x_hi;
	double y_lo;
	double y_hi;
	double[16] bs;
	Tree myTree = new Tree(0,100,100,0);
	//First two lines have some tree data
	myTree.globalMaxError = to!double(chomp(treeFile.readln())); 
	myTree.globalMinArea = to!double(chomp(treeFile.readln()));
	myTree.X_min = to!double(chomp(treeFile.readln()));
	myTree.X_max = to!double(chomp(treeFile.readln()));
	myTree.Y_min = to!double(chomp(treeFile.readln()));
	myTree.Y_max = to!double(chomp(treeFile.readln()));
	while (!treeFile.eof()){
		for (int j = 0; j != 9; j++){
			lines[j] = chomp(treeFile.readln());
		}
		idx = to!int(lines[0]);
		lefts ~= to!int(lines[1]);
		rights ~= to!int(lines[2]);
		splitID = to!char(lines[3]);
		x_lo = to!double(lines[4]);
		x_hi = to!double(lines[5]);
		y_lo = to!double(lines[6]);
		y_hi = to!double(lines[7]);
		foreach(b_i, b; split(lines[8])) bs[b_i] = to!double(b);
		if (i == 0) {
			myTree.Nodes[0] = new TreeNode(new patch(x_lo,x_hi,y_lo,y_hi));//there is a node already in there
			myTree.x_lo = x_lo;
			myTree.x_hi = x_hi;
			myTree.y_lo = y_lo;
			myTree.y_hi = y_hi;
			}
			else myTree.Nodes ~= new TreeNode(new patch(x_lo,x_hi,y_lo,y_hi));
		myTree.Nodes[i].idx = idx;
		myTree.Nodes[i].splitID = splitID;
		myTree.Nodes[i].nodePatch.bs = bs;
		i++;
	}
	//back fill the pointers now that the whole tree is constructed
	foreach(node_i,node; myTree.Nodes){
		if (node.splitID != 'N') {
			node.left = &myTree.Nodes[lefts[node_i]];
			node.right = &myTree.Nodes[rights[node_i]];
			}
		}
	writefln("Finished reading in EOS look up table from %s", filename);
	return myTree;
}


double[] linspace(double start, double stop, double n){
	double[] ys;
	for(int i = 0; i != n+1;i++){
		ys ~= (stop-start)*i/n + start;
	}
	return ys;
	}





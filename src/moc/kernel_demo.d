/*
 * file: moc_kernel_demo.d
 * location: /home/momar/Documents/Thesis/IMOC_D
 * author: Momar Hughes
 * date: 5 Apr 2015
 * description: Demo of moc_kernel.d
 */

// dmd kernel_demo.d gasdynamic.d linesearch.d kernel.d geom.d

import std.stdio, std.math;
import gasdynamic;
import kernel;
import geom; 

void main(){

    writeln("Begin moc_kernel_demo.d");
	
    writeln("Test: SetRatioOfSpecificHeats,GetRatioOfSpecificHeats...");
    writeln("Default gamma = ",GetRatioOfSpecificHeats());
    double gamma = 1.4;
    SetRatioOfSpecificHeats(gamma);
    writefln("Gamma = %g, GetRatioOfSpecificHeats = %g",
	     gamma,GetRatioOfSpecificHeats());
    writeln("...done");
    writeln();
		
    writeln("Test: SetAxiFlag,GetAxiFlag...");
    writeln("Default AxiFlag = ",GetAxiFlag());
    int axiflag = 1;
    SetAxiFlag(axiflag);
    writefln("axiflag = %d, SetAxiFlag = %d (0=OK),GetAxiFlag = %d",
	     axiflag,GetAxiFlag());
    axiflag = 2;
    SetAxiFlag(axiflag);
    writefln("axiflag = %d, SetAxiFlag = %d (0=OK),GetAxiFlag = %d",
	     axiflag,GetAxiFlag());
    writeln("...done");
    writeln();
		
    writefln("NumberOfNodes = %d,Gamma = %g,AxiFlag = %d",
	     GetNumberOfNodes(),GetRatioOfSpecificHeats(),GetAxiFlag());
    writeln("Node information...");	
    writefln("length = %d,type = %s",Node.length,typeof(Node).stringof);	
    writeln("...done");
    writeln();
		
    writeln("Test: CreateNode, GetNodeData...");	
    int id = 5;
    writefln("Create Node at id = %d, CreateNode=%d",id,CreateNode(id));
    auto result = GetNextNodeId(1);
    writefln("Search from id=1, GetNextNodeId: id = %d", result);
    result = GetNextNodeId(MAX_NODES-1);
    writefln("Search from id=max, GetNextNodeId: id = %d", result);
    GetNextNodeId(MAX_NODES-1);
    writefln("GetNodeData for Node %d = %s",id,GetNodeData(id));		
    writeln("...done");
    writeln();
	
    writeln("Test: DeleteNode...");
    DeleteNode(id);
    writefln("Delete Node %d, now Node[%d] = %s",id,id,Node[id]);	
    writeln("...done");
    writeln();
		
    writeln("Test: SetNodeData, FindNodesNear, ListNodesNear...");
    writeln("Create nodes at (1,0) (5,2) (15,34) (7,3)");
    CreateNode(-1),CreateNode(-1),CreateNode(-1),CreateNode(-1);
    writefln("There are now %d nodes",GetNumberOfNodes());
    id = GetNextNodeId(-1);
    SetNodeData(id,"X",1.0);SetNodeData(id,"Y",0.0,1);
    id = GetNextNodeId(id);
    SetNodeData(id,"X",7.0);SetNodeData(id,"Y",3.0);	
    id = GetNextNodeId(id);
    SetNodeData(id,"X",15.0);SetNodeData(id,"Y",34.0);
    id = GetNextNodeId(id);
    SetNodeData(id,"X",5.0,2);SetNodeData(id,"Y",2.0,2);
    SetNodeData(id,"CZeroUp",3.0);
    SetNodeData(id,"T",298.15);
    SetNodeData(id,"V",600.0);
    WriteNodeData(id);
    writeln("Nodes created");
    double tol = 10.0;
    int[] theids;
    int nearcount = 2;
    int nodecount = FindNodesNear(Vector3(1.0,0.0),tol,theids,nearcount);
    writefln("FindNodesNear id=%d,tol=%g,no.=%d",1,tol,nodecount);
    writeln(theids);
    writefln("ids of Nodes near are: %s",ListNodesNear(Vector3(1.0,0.0),tol,nearcount));
    writeln("...done");
    writeln();
	
    writeln("Test: SaveNodes...");
    SaveNodes("kernel_demo.txt");
    writeln("... node data saved.");
	
    writeln("Test: LoadNodes...");
    writeln("Clear all saved node data");
    foreach(int j;0 .. 10){DeleteNode(j);}
    writefln("There are now %d nodes",GetNumberOfNodes());
    writeln("Reload node data from saved file...");
    LoadNodes("kernel_demo.txt");
    writefln("There are now %d nodes",GetNumberOfNodes());
    writeln("...done");
    writeln();
	
    writeln("Finished kernel_demo.d");
}


/*
 * file: anderson_example_11_1.d
 * location: /home/momar/Documents/Thesis/IMOC_D
 * author: Momar Hughes
 * date: 12 May 2015
 * description: Mach 2.4 nozzle design, example 11.1 
 * from Anderson 1990 Modern Compressible Flow With 
 * Historic Perspective, 2nd edition, McGraw-Hill
 */

import std.stdio, std.math, std.algorithm;
import kernel;
import wall;
import unitproc;
import gasdynamic;
import geom;
import gpath;

void main()
{
    writeln("Begin example 11.1 from Anderson, Modern Compressible Flow with Historical Perspective (1982)");
    writeln();
    SetIsentropicFlag(YES); // isentropic solution
    SetFlowState(0); // work in 0th flow state
    SetAxiFlag(NO); // axisymmetric solution
    SetRatioOfSpecificHeats(1.40);
    SetGasConstant(287.0);
    double g = GetRatioOfSpecificHeats();
    double R = GetGasConstant();
    double T0 = 3000.0;
    double P0 = 7.0e6;
    //
    writeln("Create a wall along the x-axis...");
    WallFromPoints(0,"Bezier",[Vector3(0.0,1.0),Vector3(1.0,1.0)]);
    writeln("...done");
    writeln();
    //
    writeln("Create a centred expansion fan at (0.0,1.1)...");
    int[] fanids;
    int id;
    foreach(int i;0 .. 7){
        id = 100 + i+1;
        fanids ~= id;
        CreateNode(id);
        double tpmvalue = (0.375 + i*3.0)*PI/180.0;
        double machvalue = PM2(tpmvalue,g);
        double T = T0/T0_T(machvalue,g);
        double P = P0/p0_p(machvalue,g);
        SetNodeData(id,"X",0.0);
        SetNodeData(id,"Y",1.1);
        SetNodeData(id,"theta",tpmvalue);
        SetNodeData(id,"P",P);
        SetNodeData(id,"T",T);
        SetNodeData(id,"rho",P/R/T);
        SetNodeData(id,"V",machvalue*sqrt(g*R*T));
    } // end foreach
    writefln("There now exists %d nodes",GetNumberOfNodes());
    writeln("...done");
    writeln();
    //
    writeln("Compute the first wall node radiating from the fan.");
    int old_node = CMinusWallNode(0,-1,fanids[0],-1);
    writeln("...done");
    writeln();
    //
    writeln("Compute the rest of the fan radiating down to the wall.");
    int new_node;
    int last_node;
    int axis_node;
    foreach(i; 1 .. 7){
        new_node = fanids[i];
        int[] nodelist;
        MarchAlongC(nodelist,old_node,new_node,"minus","down");
        old_node = nodelist[1];
        last_node = nodelist[$-1];
        axis_node = CMinusWallNode(0,-1,last_node,-1);
    } // end foreach
    writefln("There now exists %d nodes",GetNumberOfNodes());
    writeln("...done");
    writeln();  
    double x_cone = GetNodeData(axis_node).x;
    double M_cone = GetNodeData(axis_node).Mach;
    writefln("Mach cone starts at x = %g with M = %g",x_cone,M_cone);
    writeln("Put down a number of nodes along the x-axis with constant M.");
    writeln("Work back upstream from each of these nodes.");
    double dL = 0.05;
    double Mach_angle = MachAngle(M_cone);
    double dx = dL * cos(Mach_angle);
    double dy = dL * sin(Mach_angle);
    old_node = last_node;
    int old_edge = axis_node;
    int new_edge;
    foreach(int i; 0 .. 11){
        new_edge = CreateNode(-1);
        SetNodeData(new_edge,"X",x_cone + (i+1)*dx);
        SetNodeData(new_edge,"Y",1.0+(i+1)*dy);
        double T = T0/T0_T(M_cone,g);
        double P = P0/p0_p(M_cone,g);
        SetNodeData(new_edge,"P",P);
        SetNodeData(new_edge,"T",T);
        SetNodeData(new_edge,"rho",P/R/T);
        SetNodeData(new_edge,"V",M_cone*sqrt(g*R*T));
        SetNodeData(new_edge,"Theta",0.0);
        SetNodeData(new_edge,"CPlusUp",old_edge);
        SetNodeData(old_edge,"CPlusDown",new_edge);
        int[] nodelist;
        MarchAlongC(nodelist,old_node,new_edge,"Minus","Up");
        old_node = nodelist[1];
        old_edge = new_edge;
    } // end foreach
    writefln("There now exists %d nodes",GetNumberOfNodes());
    writeln("...done");
    writeln();
    //  
    writeln("Start at the last node on the fan and step along a streamline");
    writeln("until either (1) we cross the characteristic defining the start");
    writeln("of the uniform flow region or (2) the Shepard interpolation fails.");
    writeln();
    double slope = dy / dx;
    int j = 0;
    int[] sn;
    sn ~= StepStreamNode(fanids[6], -1, 0.05);
    Vector3 stream_node = GetNodeData(sn[0]).pos;
    double y_cone = 1.0 + slope * fabs(stream_node.x - x_cone);
    while(stream_node.y > y_cone && sn[j] != -1){
        ++j;
        sn ~= StepStreamNode(sn[j-1], -1, 0.05);
        if(sn[j]==-1){break;}
        stream_node = GetNodeData(sn[j]).pos;
        y_cone = 1.0 + slope * fabs(stream_node.x - x_cone);    
    } // end while
    SaveNodes("anderson_11_1.txt");
    writeln("End of exercise. Nodes saved in anderson_11_1.txt.");
    writeln();
    foreach(i; 0 .. 139){
        // WriteNodeData(i);            // THIS LINE PRINTS ALL NODE DATA IF WANTED
    } // end foreach
} // end anderson_example_11_1

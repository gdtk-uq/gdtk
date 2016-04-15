/* paver.d
 * Generate a 2D unstructured grid of quadrilateral cells from a discretized boundary.
 *
 * Heather Muir - Undergraduate Thesis, 15/4/2016
 */

module paver;

import std.stdio;
import std.math;
import std.format;
import std.conv;
import std.algorithm: canFind, sort, reverse, remove; 



/*-----------------------------------Section 1----------------------------------*/
/*-------------------------------Class definition-------------------------------*/
/*---------------------------------- 116 Lines ---------------------------------*/

class Point
{
    size_t ID;
    double x;
    double y;
    double z = 0;

    double interior_angle;
    string node_cat = "none";
    double[2] V_R; //vector to node on right
    double[2] V_L; //vector to node on left
    Point pointL;
    Point pointR;

    size_t rowID;

    this(size_t ID, double x, double y){
	this.ID = ID;
	this.x = x;
	this.y = y;
    }

    void store_node_angle(double node_angle){
	this.interior_angle = node_angle;
    }

    void store_node_cat(string node_cat){
	if (this.node_cat == "none") {
	    this.node_cat = node_cat;
	}
    }

    void store_vectors(double[2] vr, double[2] vl){
	this.V_R = vr;
	this.V_L = vl;
    }

    void store_adjacent_points(Point pointL, Point pointR){
	this.pointL = pointL;
	this.pointR = pointR;
    }

    void store_rowID(size_t rowID){
	this.rowID = rowID;
    }
    
    string print(){
	return format("row ID: %s | xyz(%s,%s,%s) \n interior angle: %s | node cat: %s",rowID,x,y,z,interior_angle,node_cat);
    }  

    static size_t ID_count;
}

/*----------------------------------------------------------*/

class UFace 
{
    //edge in 2D

    size_t face_ID;
    size_t[2] point_IDs;
    
    this(size_t face_ID, size_t[2] point_IDs){
	this.face_ID = face_ID;
	this.point_IDs = point_IDs;
    }
    string print(){
	return format("face ID: %s | contains points: [%s, %s])", face_ID, point_IDs[0], point_IDs[1]);
    }  

    static size_t face_ID_count;
}

/*----------------------------------------------------------*/

class Cell
{
    size_t cell_ID;
    string cell_type;
    size_t[] point_IDs;
    size_t[] face_IDs;
    int[] outsigns; //order of outsigns corresponds to face_IDs

    this(size_t cell_ID){
	this.cell_ID = cell_ID;
    }

    void store_point_IDs(Point[] points){
	foreach(p;points){
	    this.point_IDs ~= p.ID;
	}
    }

    void store_face_IDs(UFace[] faces){
	foreach(f;faces){
	    this.face_IDs ~= f.face_ID;
	}
    }

    void assign_outsigns(int[] outsigns){
	foreach (o;outsigns){
	    assert(o == -1 || o == 1);
	}
	this.outsigns ~= outsigns;
    }

    void auto_cell_type(){
	int num_faces = to!int(this.face_IDs.length);
	auto cell_type = cast(CellType)num_faces;
	string cell = to!string(cell_type);
	this.cell_type = cell;
    }

    string print(){
	return format("cell ID: %s | cell type: %s | contains point IDs: %s and face IDs: %s",cell_ID,cell_type, point_IDs, face_IDs);
    }

    static size_t cell_ID_count;

}

/*---------------------------------------------------------*/

enum CellType { //**this is different to USGCell_type in usgrid.d
    none = 0,
    line = 1,
    dubline = 2,
    triangle = 3,
    quad = 4, 
    pent = 5,
    hex =6   
}


/*-----------------------------------Section 2----------------------------------*/
/*-----------------------gloally accessible data storage------------------------*/
/*----------------------------------- 6 Lines ----------------------------------*/

int ndimensions = 2;
Point[] POINT_LIST;
UFace[] FACE_LIST;
Cell[] CELL_LIST;
double initial_size;
double avr_area;

/*---------------------------------------*/

/*-----------------------------------Section 3----------------------------------*/
/*---------------------geometry & vector operation functions--------------------*/
/*----------------------------------- 88 Lines ---------------------------------*/

@nogc  double[2][2] left_rot_matrix(double theta)
{
    //arbitrary angle rotation matrix

    double[2][2] matrix = [[cos(theta), sin(theta)],
			[-sin(theta), cos(theta)]];
    return matrix;
}

double[2][2] left90matrix = left_rot_matrix(PI/2);

@nogc double[2] vector_rot(double[2] vector, double[2][2] matrix)
{
    //rotate a vector by a given rotation angle matrix

    double[2] result;
    result[0] = (vector[0]*matrix[0][0])+(vector[1]*matrix[1][0]);
    result[1] = (vector[0]*matrix[0][1])+(vector[1]*matrix[1][1]);
    return result;
}

Point midpoint(Point a, Point b)
{
    Point mid = new Point(0, 0.5*(a.x+b.x), 0.5*(a.y+b.y));
    return mid;
}
//^garbage collection here

@nogc double point_dist(Point a, Point b)
{
    //absolute distance between 2 points

    double dist = ((a.x-b.x)^^2+(a.y-b.y)^^2)^^0.5;
    return dist;
}

@nogc double[2] vector2D(Point a, Point b)
{
    //make 2D array vector from point a to point b

    double[2] vector;
    vector[0] = b.x-a.x;
    vector[1] = b.y-a.y;
    return vector;
}

@nogc int LeftorRight(Point a, Point b, Point X)
{
    /* 
       determine if point X lies to the left or the right of line AB
       will return 1 for "left", -1 for "right" or 2 for "on line" 
    */

    double[2] AB = vector2D(a,b);
    double[2] AX = vector2D(a,X);
    double determinant = AB[0]*AX[1] - AB[1]*AX[0];
    int ans;
    if (determinant == 0) {ans = 2;}
    else if (determinant < 0) {ans = -1;}
    else if (determinant > 0) {ans = 1;}

    return ans;
}

@nogc double norm(double[2] vector)
{
    //return the length(norm) of a 2D vector

    double dist =  ((vector[0])^^2+(vector[1])^^2)^^0.5;
    return dist;
}

@nogc double[2] unit_normal_vector(Point point0, Point pointR)
{
    double[2] V_f = vector2D(point0, pointR);
    double[2] V_norm = vector_rot(V_f, left90matrix)[]/norm(V_f);
    return V_norm;
}

@nogc double[2] unit_leftV(double[2] vr, double theta)
{
    double[2][2] matrix = left_rot_matrix(theta);
    double[2] vector = vector_rot(vr, matrix);
    vector[] /= norm(vector);
    return vector;
}

/*-----------------------------------Section 4----------------------------------*/
/*----------some useful geometry functions utilized in Paving method------------*/
/*---------------------------------- 175 Lines ---------------------------------*/

double squareness(Cell cell)
{
    int n = to!int(cell.point_IDs.length);
    assert (n==4);
    double area;
    double perimeter;
    Point[] points;
    for(int i; i<4; ++i){
	points ~= POINT_LIST[cell.point_IDs[i]];
    }

    area = 0.5*abs(vector2D(points[0], points[1])[0]*vector2D(points[0], points[3])[1] - vector2D(points[0], points[1])[1]*vector2D(points[0], points[3])[0]);
    area += 0.5*abs(vector2D(points[2], points[1])[0]*vector2D(points[2], points[3])[1] - vector2D(points[2], points[1])[1]*vector2D(points[2], points[3])[0]);

    perimeter = point_dist(points[0], points[1])+point_dist(points[1], points[2])+point_dist(points[2], points[3])+point_dist(points[3], points[0]);

    double ans = 16*area/(perimeter^^2);
    return ans;
}
//^(garbage collection here)

int winding_number(Point P, Point[] surface)
{
    //if wn is zero the point lies outside the closed boundary
    int wn;  
    Point a;
    Point b;
    Point Pr = new Point(0, P.x+1, P.y);
    foreach(node;surface){
	a = node;
	b = node.pointR;  
	//point at unit vector to right of P
	if(LeftorRight(P, Pr, b)==1 && LeftorRight(P,Pr,a)==-1 && LeftorRight(a,b,P)==1){ ++wn;}
	if(LeftorRight(P, Pr, b)==-1 && LeftorRight(P,Pr,a)==1 && LeftorRight(a,b,P)==-1){ --wn;}
    }
    return wn;
}
//^(garbage collection here)

double clockwise_check(Point[] surface)
{
    //If the result is positive the curve is clockwise, 
    //if it's negative the curve is counter-clockwise. 
    //The result is (twice?) the enclosed area, with a +/- convention.
    double sum = 0;
    foreach(p;surface){
	sum += (p.x-p.pointL.x)*(p.y+p.pointL.y);
    }
    double area = sum/2;

    closed_row_check(surface);
    //if(area > 0){writeln("closed row is ordered clockwise");}
    //else {writeln("closed row is ordered anti-clockwise");}

    return area;
    //**function returns the enclosed area
}
//^(garbage collection here)

double calc_area(Point[] cell_points)
{
    double sum = 0;
    int n = to!int(cell_points.length);
    for(int i; i<n-1; ++i) { 
	sum += (cell_points[i+1].x-cell_points[i].x)*(cell_points[i+1].y+cell_points[i].y);
    }
    sum += (cell_points[0].x-cell_points[n-1].x)*(cell_points[0].y+cell_points[n-1].y);
    double area = abs(sum/2);

    return area;
}
//^(garbage collection here)

void closed_row_check(Point[] surface)
{
    //exception thrown if points in surface are not all
    //holding hands with eachother.
    //note that node properties must be assigned to check this
    foreach(point;surface){
	if (point.pointL.pointR !is point){
	    throw new Exception(format("this point is not linked correctly: \n %s", point.print()));
	}
    }
}
//^(garbage collection here)

Point[4][] intersecting_boundary(Point[] surface)
{
    Point[4][] intersections;
    //note: if length == 0, no intersections found
    //node_properties must have been assigned
    foreach(point;surface){
	foreach(otherpoint;surface){
	    if (otherpoint !is point && otherpoint !is point.pointR && otherpoint !is point.pointL) {
		Point a = otherpoint;
		Point b = otherpoint.pointR;
		if (LeftorRight(point, point.pointR, a) == -1*LeftorRight(point, point.pointR, b) && LeftorRight(a, b, point) == -1*LeftorRight(a,b,point.pointR)) {
		    Point[4] points = [point, point.pointR, a, b];
		    intersections ~= points;
		}
	    }
	}
    }
    return intersections;
}
//^(garbage collection here)

Point[4][] intersecting_faces(Point[] surface)
{
    //note: specifically for checking after each new row is generated
    Point[4][] intersections;
    size_t npoints_row = surface.length;
    UFace[] new_row_faces = FACE_LIST[$-npoints_row..$];
    Point a;
    Point b;
    Point c;
    Point d;
    //note: if length == 0, no intersections found
    //node_properties must have been assigned
    foreach(face;new_row_faces){
	foreach(otherface;FACE_LIST){
	    if (otherface !is face) {
	        a = POINT_LIST[face.point_IDs[0]];
	        b = POINT_LIST[face.point_IDs[1]];
	        c = POINT_LIST[otherface.point_IDs[0]];
	        d = POINT_LIST[otherface.point_IDs[1]];
		if (LeftorRight(a,b,c) == -1*LeftorRight(a,b,d) && LeftorRight(c,d,a) == -1*LeftorRight(c,d,b)) {
		    Point[4] points = [a,b,c,d];
		    intersections ~= points;
		}
	    }
	}
    }
    return intersections;
}
//^(garbage collection here)

@nogc double average_face_length(Point[] surface)
{
    double sum = 0;
    int count;
    foreach(point;surface){
	++count;
	sum += point_dist(point, point.pointR);
    }
    double ans = sum/count;
    return ans;
}

double angle_between_vectors(double[2] a, double[2] b)
{
    double dot = a[0]*b[0]+a[1]*b[1];
    double abs = norm(a)*norm(b);
    double theta = acos(dot/abs);
    Point pa = new Point(0, 0, 0);
    Point pb = new Point(0, a[0], a[1]);
    Point X = new Point(0, b[0], b[1]);
    if (LeftorRight(pa,pb,X) != -1){
	theta = 2*PI - theta;
    }
    return theta;
}
//^(garbage collection here)

bool inside_prev_row_check(Point[] New_Row, Point[] Prev_Row){
    bool inside = 1;
    int wn;
    foreach(p; New_Row){
	wn = winding_number(p, Prev_Row);
	if(wn==0){
	    inside = 0;
	    break;
	}
    }
    return inside;
}
//^(garbage collection here)

/*---------------------------------------------------------------------------*/

/*-----------------------------------Section 5----------------------------------*/
/*----------------------------MAJOR PAVING FUNCTIONS----------------------------*/
/*---------------------------------- 1018 Lines --------------------------------*/

Point make2DPoint(double x, double y)
{
    size_t ID = Point.ID_count;
    ++Point.ID_count;

    Point point = new Point(ID, x, y);
    POINT_LIST ~= point;

    return point;
}

UFace makeFace(size_t ID_0, size_t ID_1)
{
    /* in 2D this is a line from point ID_0
       to point ID_1 */

    size_t ID = UFace.face_ID_count;
    ++UFace.face_ID_count;

    UFace face = new UFace(ID, [ID_0, ID_1]);
    FACE_LIST ~= face;

    return face;
}

Cell makeCell(Point[] points)
{
    //**NOTE** points must be passed in anticlockwise order!
    // or outsigns will be calculated incorrectly. 

    size_t ID = Cell.cell_ID_count;
    ++Cell.cell_ID_count;
    UFace[] faces;
    int[] outsigns;

    for(uint i; i<points.length; ++i){
	auto a = i;
	uint b;
	if (a == points.length-1) {b = 0;}
	else {b = a+1;}
	UFace found1 = findFace([points[a].ID, points[b].ID], 1);
	UFace found2 = findFace([points[b].ID, points[a].ID], 1);
	if (found1 is null && found2 is null){
	    faces ~= makeFace(points[a].ID, points[b].ID);
	    outsigns ~= -1;
	} else if (found1 is null) {
	    faces ~= found2;
	    outsigns ~= 1;
	} else {
	    faces ~= found1;
	    outsigns ~= -1;
	}
    }

    Cell cell = new Cell(ID);
    cell.store_point_IDs(points);
    cell.store_face_IDs(faces);
    cell.assign_outsigns(outsigns);
    cell.auto_cell_type();
    //could adapt this class to make these all part of the constructor
    CELL_LIST ~= cell;

    return cell;
}

UFace findFace(size_t[2] pointIDs, bool directional)
{
    /*
      From 2 point IDs find the face in the global FACE_LIST.
      If directional == 0 test both orders.
      If face does not exist, function returns null.
    */

    UFace found;
    size_t[2] reverse_case;
    if (directional == 0){
	reverse_case = [pointIDs[1], pointIDs[0]];
    }
    foreach(face;FACE_LIST){
	if(face.point_IDs == pointIDs || face.point_IDs == reverse_case){
	    found = face;
	}
    }
    return found;
}

void assign_node_properties(Point[] surface)
{
    /* 
       'surface' is the list of points which make up the 
       current paving boundary.
       
       For each node in the surface, function calculates and assigns:
       -interior angle --> node category
       -adjacent points (PointL & PointR)
       -vectors to attached points (V_L & V_R)
    */

    size_t n = surface.length;

    for(size_t i; i<n; ++i){

	size_t nodeID = i;
	Point point0 = surface[nodeID];  
	Point pointR; //must be created outside of conditional statements
	Point pointL;

	if(nodeID==n-1){
	    pointL = surface[nodeID-1];
	    pointR = surface[0];
	} else if(nodeID==0){
	    pointL = surface[n-1];
	    pointR = surface[nodeID+1];
	} else {
	    pointL = surface[nodeID-1];
	    pointR = surface[nodeID+1];
	}

	/*------------Diagram--------------
                      
	       R----<----0 
	                 |  (exterior)
	     (interior)  |  
	                 ^  <- (indicates anti-clockwise
	                 |             order)
	                 |
	                 L
	  --------------------------------*/

	//good old cosine rule:
	double c = point_dist(pointL, pointR);
	double a = point_dist(pointL, point0);
	double b = point_dist(pointR, point0);
	double r = (a^^2+b^^2-c^^2)/(2*a*b);
	assert(r<1.01 && r>-1.01); //checking calculated r is plausible
	if (r >1.0) r=1.0;
	if (r < -1.0) r=-1.0;
	double angle = acos(r);

	if (angle > 3.141 && angle < 3.142){
	    //angle is basically pi it doesn't need adjusting
	} else {
	    // determine if the interior or exterior angle was calculated
	    // ie. cosine rule will just calculated smallest angle
	    // see if point0 is on left or right of line L-->R
	    // point0 should be to the right for interior angle calc

	    int ans = LeftorRight(pointL, pointR, point0);
	    // function returns 1 for "left", -1 for "right" or 0 for "on line"

	    if (ans == 1) { //exterior angle was calculated - fix:
		angle = 2*PI - angle;
	    }
	}

	//store the node information:
	point0.store_node_angle(angle);
	string node_cat = node_category(angle);
	point0.store_node_cat(node_cat);
	double[2] V_R = vector2D(point0, pointR);
	double[2] V_L = vector2D(point0, pointL);
	point0.store_vectors(V_R, V_L);
	point0.store_adjacent_points(pointL, pointR);

    } 
}

void adjust_row_ends(Point[] surface, double initial_size,  bool boundary)
{
    /*
      prupose 1:
      If neighbouring row ends appear, the smaller angle of the two
      nodes wins the battle for "row end" status.
      The label "with row end" is assigned to the node on the left
      of row ends, which will close off the cell, unless another row end 
      exists 2 nodes to the left (causing consecutive end/with-end/end) 
      when the node attached to the right side is granted "with row end"
      status. A few other special cases which cause issues are handled 
      as well.
      purpose 2:
      Either if too few row ends exist to or if cell size is still 
      reducing too severly from the initial_size, additional row ends 
      are introduced to manage this cell size reduction.
      This is disabled for the first row generated from the boundary,
      ie. when bool boundary == 1.
    */

    bool conditions(Point node){
	bool ans;
	if(node is null){return ans;}
	if(node.pointL.node_cat != "row end" && 
	   node.pointR.node_cat != "row end" && 
	   node.pointL.pointL.node_cat != "row end" && 
	   node.pointR.pointR.node_cat != "row end" && 
	   node.pointL.node_cat != "with row end" && 
	   node.pointR.node_cat != " with row end" &&
	   node.pointL.pointL.node_cat != "with row end" && 
	   node.pointR.pointR.node_cat != " with row end"){
	    ans = 1;
	} else{ ans = 0;}
	return ans;
    }

    //interesting to see whether this section helps:
    foreach(node;surface){
	size_t[] point_IDs;
	size_t node_count;
	foreach(face; FACE_LIST){
	    point_IDs = face.point_IDs;
	    if (point_IDs.canFind(node.ID)){
		++node_count;
	    }
	}
	if (node_count >= 4 && conditions(node)){
	    node.node_cat = "row end";
	}
    }
    //----------------------------------------------
    
    if (boundary == 0){
	Point[] row_end_list;
	double smallest_angle = 2*PI;
	Point smallest_angle_node;

	foreach (node;surface){
	    if (node.node_cat == "row end") {
		row_end_list ~= node;
	    } else if (node.interior_angle < smallest_angle){
		smallest_angle = node.interior_angle;
		smallest_angle_node = node;
	    }
	}
	
	//writefln("there were %s row ends in this row", row_end_list.length);
	if (row_end_list.length == 0){
	    smallest_angle_node.node_cat = "row end";
	    row_end_list ~= smallest_angle_node;
	}
	if (row_end_list.length == 1){
	    size_t move = surface.length/2;
	    Point node = smallest_angle_node;
	    for(size_t i; i<move; ++i){
		node = node.pointR;
	    }
	    if (node.interior_angle < PI){
		node.node_cat = "row end";
		row_end_list ~= node;
	    } else {//throw new Exception(text("the node we wanted to enforce as a matching row end had too big of an angle... might need to deal with this?"));
	    }
	}
        
	double avr_length = average_face_length(surface);
	size_t n = surface.length;
	//writeln("comparing: ", avr_length*n/(n-2), " with ", prev_avr_length);
	if (avr_length*n/(n-2) < initial_size){
	    double smallest = 2*initial_size;
	    double biggest = 2*avr_length;
	    double size;
	    Point smallest_node;
	    Point biggest_node;
	    foreach(node; surface){
		size = point_dist(node, node.pointR) + point_dist(node, node.pointL);
		if (size < smallest && node.node_cat != "row end"){
		    smallest = size; smallest_node = node;
		} else if (size > biggest){
		    biggest = size; biggest_node = node;
		}
	    }
	    
	    if(conditions(smallest_node)){ 
		writeln("enforcing row end to help manage cell size reduction"); 
		smallest_node.node_cat = "row end";
	    } else if(conditions(smallest_node.pointL)){
		writeln("enforcing row end to help manage cell size reduction");
		smallest_node.pointL.node_cat = "row end";
	    } else if(conditions(smallest_node.pointR)){
		writeln("enforcing row end to help manage cell size reduction");
		smallest_node.pointR.node_cat = "row end";
	    } 
	    if (biggest > 4*avr_length){ 
		//writeln("big detected"); biggest_node.node_cat = "row corner";
		// currently only the small case is dealt with.
		// the large cell case is better dealt with in 'fix_big_cells'
	    }
	}
    }  

    foreach(node; surface){
	if (node.node_cat == "row end" && node.pointR.node_cat == "row end"){
	    if (node.interior_angle < node.pointR.interior_angle){
		node.pointR.node_cat = "row side";
	    }
	    else if(node.interior_angle >= node.pointR.interior_angle){
		node.node_cat = "row side";
	    }
	}
	if(node.pointL.pointL.node_cat == "row end" &&
	   node.pointR.pointR.node_cat == "row end"){
	    node.node_cat = "row side";
	}
	if (node.node_cat == "row end") {
	    if (node.pointL.pointL.node_cat != "row end"){
		node.pointL.node_cat = "with row end";
	    } else {
		node.pointR.node_cat = "with row end";
		node.pointL.node_cat = "row side";
	    }
	}
    }
}


string node_category(double angle)
{
    /*
      categorisation of nodes based on interior angle
      **NOTE: not dealing with ambiguous cases yet
    */

    double[6] tols = 0.35; //all tols 20 degrees for now

    string node_cat;
    if (angle <= PI/2+tols[0]) {node_cat = "row end";}
    else if ((PI/2 + tols[0] < angle) && (angle <= PI-tols[1])) {
	//node_cat = "ambiguous row end/side"; **coming soon 
	node_cat = "row side";
    }
    else if ((PI - tols[1] < angle) && (angle <= PI + tols[2])) {
	node_cat = "row side";
    }
    else if ((PI + tols[2] < angle) && (angle <= 3*PI/2 - tols[3])) {
	//node_cat = "ambiguous row side/corner"; **
	node_cat = "row side";
    }
    else if ((3*PI/2 - tols[3] < angle ) && (angle <= 3*PI/2 + tols[4])) {
	node_cat = "row corner";
    }
    else if ((3*PI/2 + tols[4] < angle) && (angle <= 2*PI - tols[5])) {
	//node_cat = "ambiguous row corner/reversal"; **
	node_cat = "row corner";
    }
    else if (angle > 2*PI) {
	throw new Error(text("calculated angle > 2*PI... error"));
    }
    else {
	node_cat = "row reversal";
    }

    return node_cat;
}


void node_gen(Point[] surface, size_t i, ref Point[] New_Nodes)
{

    /*
      Generation of new nodes in the new row, based on current row
      node categorisation and relative position.
    */

    Point node = surface[i];
    
    if (node.node_cat == "row end" || node.node_cat == "with row end"){
	// no new nodes are generated
    }
    if (node.node_cat == "row side"){
	// generates 1 new node
	// makes 1 new cell

	/*----------------diagram:------------------
	
	            ............o  <-- (new node)
                    .           |   
		    .   new     | 
		    .   cell    |
		    .           |
       (pointL) >   o-----------o-----------o  < (pointR)
	                        ^ 
		          (current node)

	   -----------------------------------------*/

	double V_length = 0.5*(norm(node.V_L)+norm(node.V_R)/sin(node.interior_angle/2));
	double[2] node_V = [node.x, node.y];
	double[2] unit_V = unit_leftV(node.V_R, node.interior_angle/2);
	double[2] node_pos = node_V[]+(V_length*unit_V[])[];
	New_Nodes ~= make2DPoint(node_pos[0], node_pos[1]);

	if(node.pointL.node_cat == "row end") {
	    makeCell([New_Nodes[$-1], node.pointL.pointL, node.pointL, node]);
	    New_Nodes[$-1].node_cat = "row end";
	    if (New_Nodes.length >= 2) {
		makeCell([New_Nodes[$-1], New_Nodes[$-2], node.pointL.pointL.pointL, node.pointL.pointL]);
	    } 
	} else if(node.pointL.node_cat == "with row end"){
	    if (New_Nodes.length >=2 ){
		makeCell([New_Nodes[$-2], node.pointL.pointL.pointL, node.pointL.pointL, node.pointL]);
		New_Nodes[$-2].node_cat = "row end";
		makeCell([New_Nodes[$-1], New_Nodes[$-2], node.pointL, node]);
	    }	    
	} else if(New_Nodes.length >=2) {
	    makeCell([New_Nodes[$-1], New_Nodes[$-2], node.pointL, node]);
	}
		       
    }
    if (node.node_cat == "row corner") {
	// generates 3 new nodes
	// creates 2 new cells

	/*----------------diagram:------------------
	
	                   (new node1)
	            ............o...........o (new node2)
                    .           |           .
		    .   new     |    new    .
		    .  cell1    |   cell2   .
		    .           |           .
       (pointL) >   o-----------o-----------o (new node3)
	                      / |           
		       current	|        
		        node	|      
				|           
				o
				^
			     (PointR)

	   -----------------------------------------*/

	double V_length = 0.5*(norm(node.V_L)+norm(node.V_R)/sin(node.interior_angle/3));
	double[2] node_V = [node.x, node.y];
	double[2] unit_V = unit_leftV(node.V_R, 2*node.interior_angle/3);
	double[2] node_pos = node_V[]+(V_length*unit_V[])[];
        New_Nodes ~= make2DPoint(node_pos[0], node_pos[1]); //new node 1
	double V2_length = SQRT2*V_length;
	unit_V = unit_leftV(node.V_R, node.interior_angle/2);
	node_pos = node_V[]+(V2_length*unit_V[])[];
        New_Nodes ~= make2DPoint(node_pos[0], node_pos[1]); //new node 2
	unit_V = unit_leftV(node.V_R, node.interior_angle/3);
	node_pos = node_V[]+(V_length*unit_V[])[];
        New_Nodes ~= make2DPoint(node_pos[0], node_pos[1]); //new node 3

	if(node.pointL.node_cat == "row end"){
	    makeCell([New_Nodes[$-3], node.pointL.pointL, node.pointL, node]);
	} else {
	    //New_Nodes list should not be empty since start node is a row side 	    //ie. New_Nodes[$-4] should be accessible
	    makeCell([New_Nodes[$-3], New_Nodes[$-4], node.pointL, node]);
	}
	makeCell([New_Nodes[$-1], New_Nodes[$-2], New_Nodes[$-3], node]);
    }
    if (node.node_cat == "row reversal") {
	// generates 5 new nodes
	// creates 3 new cells

	/*----------------diagram:------------------
	
	                   (new node3)
      (new node2) > o...........o...........o (new node4)
                    .           |           .
		    .   new     |    new    .
		    .  cell2    |   cell3   .
		    .           |           .
     (new node1) >  o-----------o-----------o (new node5)
	           /           / \
		  /   new     /   \
		 /   cell1   /     \
		/	    /       \    
		o----------o	     o
			   ^	     ^
		       (PointL)   (PointR)

	   -----------------------------------------*/

	double V_length = 0.5*(norm(node.V_L)+norm(node.V_R)/sin(node.interior_angle/4));
	double[2] node_V = [node.x, node.y];
	double[2] unit_V = unit_leftV(node.V_R, 3*node.interior_angle/4);
	double[2] node_pos = node_V[]+(V_length*unit_V[])[];
        New_Nodes ~= make2DPoint(node_pos[0], node_pos[1]); //new node 1
	double V2_length = SQRT2*V_length;
	unit_V = unit_leftV(node.V_R, 5*node.interior_angle/8);
	node_pos = node_V[]+(V2_length*unit_V[])[];
        New_Nodes ~= make2DPoint(node_pos[0], node_pos[1]); //new node 2
	unit_V = unit_leftV(node.V_R, node.interior_angle/2);
	node_pos = node_V[]+(V_length*unit_V[])[];
        New_Nodes ~= make2DPoint(node_pos[0], node_pos[1]); //new node 3
	unit_V = unit_leftV(node.V_R, 3*node.interior_angle/8);
	node_pos = node_V[]+(V2_length*unit_V[])[];
        New_Nodes ~= make2DPoint(node_pos[0], node_pos[1]); //new node 4
	unit_V = unit_leftV(node.V_R, node.interior_angle/4);
	node_pos = node_V[]+(V_length*unit_V[])[];
        New_Nodes ~= make2DPoint(node_pos[0], node_pos[1]); //new node 5

	if(node.pointL.node_cat == "row end"){
	    makeCell([New_Nodes[$-5], node.pointL.pointL, node.pointL, node]);
	    //row end right beside row reversal
	}
	else {
	    //New_Nodes list should not be empty since start node is a row side 
	    //ie. New_Nodes[$-6] should be accessible
	    makeCell([New_Nodes[$-5], New_Nodes[$-6], node.pointL, node]);
	}
	makeCell([New_Nodes[$-3], New_Nodes[$-4], New_Nodes[$-5], node]);
	makeCell([New_Nodes[$-1], New_Nodes[$-2], New_Nodes[$-3], node]);
    }
}


Point[] generate_new_row(Point[] surface)
{
    /*
      as the name suggests, generates the new row of nodes 
      (calling on node_gen function)
    */

    Point[] New_Nodes;
    size_t n = surface.length;

    // the start node will be the second consecutive side node:
    /* note that a continuation condition is that the row contains
       at least 2 consecutive side nodes, but an exception is 
       included to ensure this has been correctly enforced */
    Point start_node;
    foreach(node;surface){
	if(node.node_cat == "row side" && node.pointL.node_cat == "row side") {
	    start_node = node;
	    break;
	}
    }
    if(start_node is null){throw new Exception("there are not 2 consecutive side nodes in this row");}

    //construct order from start node
    for(size_t i = start_node.rowID; i < n; ++i){
	node_gen(surface, i, New_Nodes);
    }
    for(size_t i = 0; i < start_node.rowID; ++i){
	node_gen(surface, i, New_Nodes);
    }

    // make final cell (known to be this case becasue of 
    // second side node start condition)
    makeCell([surface[start_node.rowID], New_Nodes[0], New_Nodes[$-1], start_node.pointL]);

    //assigning row IDs to the nodes in the new row
    size_t rowID;
    foreach(node;New_Nodes){
	node.store_rowID(rowID);
	++rowID;
    }

    return New_Nodes;

}

double summed_squareness(Point node)
{
    /*
      utilised by the smoothing function
    */

    int cell_count;
    double summed_squareness;
    foreach(cell; CELL_LIST) {
	size_t[] ID_list;
	ID_list ~= cell.point_IDs;
	if(ID_list.canFind(node.ID)) {
	    ++cell_count;
	    summed_squareness += squareness(cell);
	}
    }
    summed_squareness /= cell_count;
    return summed_squareness;
}
    


double smoothing(Point[] surface, Point[] Prev_Row, bool ghost)
{
    /* 
       Implementing Lapacian smoothing: this method moves each node
       one by one based on the average of the vectors constructed to
       every other node it attaches to via a face.
       Improved by encorporating a function which also influences the 
       movement of the node based on improving the summed_squareness 
       of all of the cells neighbouring the node.    
    */ 

    //surface ~= surface[0];

    double move_length;
    double biggest_move = 0;
    Point attached_point;
    double x;
    double y;
    double prev_squareness;
    double squareness;

    foreach(node;surface) {
	if(node.node_cat != "row corner" && node.node_cat != "row reversal"){
	    double[2][] vectors;
	    foreach(face; FACE_LIST) {
		size_t[] ID_list;
		ID_list ~= face.point_IDs;
		if(ID_list.canFind(node.ID)) {
		    if (node.ID == face.point_IDs[0]) {
			attached_point = POINT_LIST[face.point_IDs[1]];
		    } else {
			attached_point = POINT_LIST[face.point_IDs[0]];
		    }
		    if (ghost == 1){
			if(!Prev_Row.canFind(attached_point)){
			    vectors ~= [attached_point.x-node.x, attached_point.y-node.y];
			}
		    } else if (ghost == 0){
			vectors ~= [attached_point.x-node.x, attached_point.y-node.y];
		    }
		}
	    }

	    if (vectors.length == 0) {
		//writeln("apparently this node has no friends! \n", node.print());
		//throw new Error(text("apparently this node has no friends!"));
	    }
	
	    if(vectors.length != 0){
		double[2] movement = vectors[0];
		for (uint i = 1; i < vectors.length; ++i){
		    movement[] += vectors[i][];
		}
		movement[] /= 2*to!int(vectors.length); 
		move_length = norm(movement);
		if (move_length > biggest_move) { biggest_move = move_length;}

		prev_squareness = summed_squareness(node);
		x = node.x;
		y = node.y;
		node.x += movement[0];
		node.y += movement[1];
		squareness = summed_squareness(node);
		if(squareness < prev_squareness){
		    node.x = x;
		    node.y = y;
		}
	    }
	}
    }

    return biggest_move;

}

Point[] seam(Point start, ref Point[] New_Row, ref size_t[] adjust_IDs, ref size_t original_nodes_left, bool closing)
{
    /*
      constructs points for 'seaming cells' to close regions with a 
      continued acute boundary angle, or when completing the closing 
      seam (bool closing == 1).
      This process involves lots of complicated intracies - which would
      take to much space to explain here.
    */ 

    //size_t original_nodes_left = New_Row.length;
    Point[] New_Nodes;
    Point[] empty;
    Point[] inner_nodes;
    //special case - only 1 cell remianing to close
    if (closing == 1 && New_Row.length == 4){
	makeCell([start, start.pointR, start.pointR.pointR, start.pointL]);
	writeln("CLOSED (special)");
	return New_Nodes;
    }

    double midx = (start.pointL.x+start.pointR.x)/2;
    double midy = (start.pointL.y+start.pointR.y)/2;
    New_Nodes ~= make2DPoint(midx, midy);
    original_nodes_left -= 3;
    makeCell([New_Nodes[0], start.pointL, start, start.pointR]);

    Point nextpointL = start.pointL;
    Point nextpointR = start.pointR;

    double[2] a = vector2D(nextpointL, nextpointL.pointL);
    double[2] b = vector2D(nextpointR, nextpointR.pointR);
    double theta = angle_between_vectors(a,b);
    int loop;

    bool conditions(double theta, Point nextpointL, 
		    Point nextpointR, size_t original_nodes_left){
	bool conditions;
	if(closing == 0){ //not doing the closing seam
	    if(theta < PI/2 && 
	  nextpointL.interior_angle > 2.4 && 
	  nextpointL.interior_angle < 3.9 &&
	  nextpointR.interior_angle > 2.4 &&
	  nextpointR.interior_angle < 3.9 && 
	  original_nodes_left > 2 &&
	  !adjust_IDs.canFind(nextpointL.rowID) &&
	  !adjust_IDs.canFind(nextpointR.rowID)) { conditions = 1;}
	    else{
		writefln("there are %s remaining nodes \n (not closing seam)", original_nodes_left);
		conditions = 0;}
	} else if(closing == 1){
	    if(original_nodes_left > 2){conditions = 1;}
	    else{
		writeln("less than 2 remaining nodes");
		writefln("there are %s remaining nodes", original_nodes_left);
		conditions = 0;}
	}
	return conditions;
    }
	     
    int ind = 2; //reference index for joining (complicated)
    int nodes = 1;
    int Flag = 1;
    while(conditions(theta, nextpointL, 
		     nextpointR, original_nodes_left)){
	++loop;
	double avr_length = (norm(a)+norm(b))/2;
	if (point_dist(nextpointL.pointL, nextpointR.pointR) < 3*initial_size){
	    if (Flag == 1){
		midx = (nextpointL.pointL.x+nextpointR.pointR.x)/2;
		midy = (nextpointL.pointL.y+nextpointR.pointR.y)/2;
		New_Nodes ~= make2DPoint(midx, midy);
		original_nodes_left -=2;
		makeCell([New_Nodes[$-1], nextpointL.pointL, nextpointL, New_Nodes[$-ind]]);
		makeCell([New_Nodes[$-1], New_Nodes[$-ind], nextpointR, nextpointR.pointR]);
		ind = 2;
		nodes = 1;
	    } else if (Flag == 2){
		ind = 4;
		makeCell([New_Nodes[$-3], nextpointL.pointL, nextpointL, New_Nodes[$-2]]);
		makeCell([New_Nodes[$-3], New_Nodes[$-1], nextpointR, nextpointR.pointR]);
		original_nodes_left -=2;
		nodes = 0;
	    }
	    Flag = 1;
	} else if(closing==1 || point_dist(nextpointL.pointL, nextpointR.pointR)<5*initial_size){
	    midx = (nextpointL.pointL.x+nextpointR.pointR.x)/2;
	    midy = (nextpointL.pointL.y+nextpointR.pointR.y)/2;
	    New_Nodes ~= make2DPoint(midx, midy);
	    midx = (nextpointL.pointL.x+New_Nodes[$-1].x)/2;
	    midy = (nextpointL.pointL.y+New_Nodes[$-1].y)/2;
	    New_Nodes ~= make2DPoint(midx, midy);
	    midx = (New_Nodes[$-2].x+nextpointR.pointR.x)/2;
	    midy = (New_Nodes[$-2].y+nextpointR.pointR.y)/2;
	    New_Nodes ~= make2DPoint(midx, midy);
	    original_nodes_left -=2;
	    if (Flag == 1){
		if (nodes==0){
		    ind = 6;
		} else {ind=4;}
		makeCell([New_Nodes[$-3], New_Nodes[$-2], New_Nodes[$-ind], New_Nodes[$-1]]);
		makeCell([New_Nodes[$-2], nextpointL.pointL, nextpointL, New_Nodes[$-ind]]);
		makeCell([New_Nodes[$-1], New_Nodes[$-ind], nextpointR, nextpointR.pointR, ]);
		ind=4;
	    } else if (Flag == 2){
		makeCell([New_Nodes[$-2], nextpointL.pointL, nextpointL, New_Nodes[$-5]]);
		makeCell([New_Nodes[$-3], New_Nodes[$-2], New_Nodes[$-5], New_Nodes[$-6]]);
		makeCell([New_Nodes[$-1], New_Nodes[$-3], New_Nodes[$-6], New_Nodes[$-4]]);
		makeCell([New_Nodes[$-1], New_Nodes[$-4], nextpointR, nextpointR.pointR, ]);
	    }

	    Flag = 2;
	    nodes = 3;
	} else { writeln("more than 5 cell spaces between joining nodes");
	    if (closing == 0){break;}
	    else { throw new Error(text("cannot terminate at >5 node seam for closing procedure - try recursion for closing"));}
	}

	nextpointL = nextpointL.pointL;
	nextpointR = nextpointR.pointR;
	a = vector2D(nextpointL, nextpointL.pointL);
	b = vector2D(nextpointR, nextpointR.pointR);
	theta = angle_between_vectors(a,b);
	//writeln("theta: ", theta);
    } 
    if (nodes == 1){
	inner_nodes ~= New_Nodes[$-1];
    } else if (nodes == 3){
	inner_nodes = [New_Nodes[$-2]]~[New_Nodes[$-3]]~[New_Nodes[$-1]];
    }
    
    if(original_nodes_left == 1 && closing == 1){
	writeln("closing final cells");
	if(Flag == 1){
	    makeCell([nextpointR.pointR, nextpointL, New_Nodes[$-(ind-1)], nextpointR]);
	    writeln("CLOSED1");
	} else if(Flag == 2){
	    makeCell([nextpointR.pointR, nextpointL, New_Nodes[$-2], New_Nodes[$-3]]);
	    makeCell([nextpointR.pointR, New_Nodes[$-3], New_Nodes[$-1], nextpointR]);
	    writeln("CLOSED2");
	}
    }
        
    adjust_IDs ~= [nextpointL.rowID, nextpointR.rowID];
    return inner_nodes;
}

void closing_seam(Point[] Final_Row)
{
    /*
      some necessary adjustments to be able to use the 'seam' 
      function to complete the closing procedure for the mesh
    */

    size_t[] adjust_IDs; // don't actually use this, but I think I need it 
                         // in this scope to call the function here. 
    assert(Final_Row.length%2 == 0, "final boundary not even number of nodes - this is not ok... and shouldn't have happened");
    Point start = Final_Row[0];
    foreach(node; Final_Row){
	if(node.interior_angle < start.interior_angle){
	    start = node;
	}
    }
    size_t nodes_left = Final_Row.length;
    seam(start, Final_Row, adjust_IDs, nodes_left, 1);
}


void fix_IDs(UFace[] face_list, Cell[] cell_list)
{
    /* 
       the row_ID numbering can't get diturbed after the seaming function
       generates nodes which are not part of a row. This function corrects
       for that.
    */

    size_t index = 0;
    foreach(face; face_list){
	face.face_ID = index;
	++index;
    }
    index = 0;
    foreach(cell;cell_list){
	cell.cell_ID = index;
	++index;
    }
}

void seam_big_cells(in size_t npoints_boundary, in double initial_size, bool fix_meganodes)
{
    /*
      THERE ARE STILL SOME BUGS IN HERE
      So far this has been determined to be the best way to reduce the size 
      of overly large cells in the constructed mesh.
      when fix_meganodes is true the function also eliminates nodes with 
      6 or more attached faces. However, this often creates other meganodes
      and the problem is just propagated. 
    */

    Point worst;
    UFace[] worst_faces;
    double worst_avr = initial_size;
    writeln("initial: ", initial_size);
    foreach(node; POINT_LIST[npoints_boundary..$]){
	UFace[] faces;
	double avr = 0;
	size_t[] point_IDs;
	foreach(face; FACE_LIST){
	    point_IDs = face.point_IDs;
	    if (point_IDs.canFind(node.ID)){
		faces ~= face;
	    }
	}
	if(fix_meganodes == 1){
	    if (faces.length >= 6){
		worst = node; worst_faces = faces;
		writeln("meganode found!");
		break;
	    } 
	}
	foreach(face; faces){
	    avr += point_dist(POINT_LIST[face.point_IDs[0]], POINT_LIST[face.point_IDs[1]])/faces.length;
	}
	//if(node.ID%10 == 0){writeln(avr); writeln(faces.length);}
	if (avr > worst_avr){
	    worst_avr = avr; worst = node; worst_faces = faces;
	}
    }
    if (worst is null){writeln("did not find any big nodes");}
    writefln("coords of worst: (%s, %s)", worst.x, worst.y);
    if(worst_faces.length < 6){fix_meganodes = 0;}

    Cell[] worst_cells;
    size_t[] worst_cell_IDs;
    Point[] worst_points;
    size_t[] worst_IDs;
    foreach(cell; CELL_LIST[npoints_boundary..$]){
	if (cell.point_IDs.canFind(worst.ID) && 
	    !worst_cell_IDs.canFind(cell.cell_ID)){
	    worst_cells ~= cell;
	    worst_cell_IDs ~= cell.cell_ID;
	    foreach(id; cell.point_IDs){
		if (id != worst.ID && !worst_IDs.canFind(id)){
		    worst_points ~= POINT_LIST[id];
		    worst_IDs ~= id;
		}
	    }
	}
    }
    Point[] ordered_points;
    ordered_points ~= worst_points[0];
    writeln(worst_points.length);

    for(size_t i; i<worst_points.length-1; ++i){
	writeln(i+1," : ", ordered_points.length);
        if(ordered_points.length != i+1){
	    throw new Exception(text("irregular geometry: cannot make sense of an ordered path"));
	}
        Point[] possible_points;
	UFace face;
	foreach(point;worst_points){
	    face = findFace([ordered_points[i].ID, point.ID], 0);
	    if (face !is null){
		possible_points ~= point;
		writeln("found face");
	    }
	}
	foreach(point;possible_points){
	    if (point !is worst && LeftorRight(worst, ordered_points[i], point)==1){
		ordered_points ~= point;
		writeln("next_node_found");
		break;
	    } //closes: test node meets condition to assign to ordered list
	} //closes: finds test node
    }//closes: do for each point in new boundary

    assign_node_properties(ordered_points);

    //removing worst point, faces and cells from global lists
    size_t[] face_IDs;
    foreach(face; worst_faces){
	face_IDs ~= face.face_ID;
    }
    sort!("a>b")(face_IDs);
    foreach(id; face_IDs){
	FACE_LIST = remove(FACE_LIST, id);
    }
    UFace.face_ID_count -= face_IDs.length;

    sort!("a>b")(worst_cell_IDs);
    foreach(id; worst_cell_IDs){
	CELL_LIST = remove(CELL_LIST, id);
    }
    Cell.cell_ID_count -= worst_cell_IDs.length;

    fix_IDs(FACE_LIST, CELL_LIST);

    //re-seam the gap:
    writeln("commencing closing seam...");
    closing_seam(ordered_points);

    Point[] inside_nodes = POINT_LIST[npoints_boundary..$];
    double biggest_move = smoothing(inside_nodes, inside_nodes, 0);
    double move = biggest_move;
    for(int i; i<3; ++i){
	move = smoothing(inside_nodes, inside_nodes, 0);
    }

    writeln("ratio: ", worst_avr/initial_size); 
    
    if(fix_meganodes == 1){
	seam_big_cells(npoints_boundary, initial_size, 1);
    } //RECURSIVE FUNCTION

    else if(worst_avr > 2*initial_size){ //tolerance can be played with
	writeln("AGAIN");
	seam_big_cells(npoints_boundary, initial_size, 0);
    } //RECURSIVE FUNCTION
        
}

double biggest_cell(Cell[] cells, double initial_size)
{
    /*
      used to determine if big cell size reduction is necessary
    */

    double biggest_area = initial_size^^2;
    double current_area;
    foreach(cell;cells){
	current_area = calc_area([POINT_LIST[cell.point_IDs[0]],
				  POINT_LIST[cell.point_IDs[1]],
				  POINT_LIST[cell.point_IDs[2]],
				  POINT_LIST[cell.point_IDs[3]]]);
	if (current_area > biggest_area){
	    biggest_area = current_area;
	}
    }
    return biggest_area;
}

void even_out_size(Point[] surface)
{
    /* function moves points in current paving row to even out
       variances in size */
    size_t n = surface.length;
    Point start = surface[0];
    double biggest = norm(surface[0].V_L) + norm(surface[0].V_R);
    double size;
    foreach (node; surface){
	size = norm(node.V_L) + norm(node.V_R);
	if (size > biggest){
	    biggest = size;
	    start = node;
	}
    }

    Point currentL = start;
    Point currentR = start;
    for (int i; i<n/2;++i){
	currentL.pointL.x = (midpoint(currentL, currentL.pointL.pointL).x+
			     currentL.pointL.x)/2;
	currentL.pointL.y = (midpoint(currentL, currentL.pointL.pointL).y+
			     currentL.pointL.y)/2;
	currentR.pointR.x = (midpoint(currentR, currentR.pointR.pointR).x+
			     currentR.pointR.x)/2;
	currentR.pointR.y = (midpoint(currentR, currentR.pointR.pointR).y+
			     currentR.pointR.y)/2;
	currentL = currentL.pointL;
	currentR = currentR.pointR;
    }
}

void point_ID_check(){
    size_t npoints = POINT_LIST.length;
    for (size_t i; i<npoints; ++i) {
	if (POINT_LIST[i].ID != i) {
	    throw new Error(text("Point ID does not match position in POINT_LIST"));
	}
    }
}

/*-----------------------------------Section 6----------------------------------*/
/*----------------------------MAJOR PAVING PROCEDURES---------------------------*/
/*---------------------------------- 195 Lines ---------------------------------*/

void pave_rows_until_intersection(ref Point[] New_Row, ref Point[] Prev_Row, in size_t npoints_boundary, ref bool seam_next)
{
    double biggest_move;
    double move;
    Point[] inside_nodes;

    size_t npoints_row;
    size_t npoints_new_row;

    Point[4][] intersections;
    int loop_counter;
    double remaining_area;
    bool side_node_Flag = 0;
    bool inside_prev;

    Point[] point_list_backup;
    size_t point_count_backup;
    UFace[] face_list_backup;
    size_t face_count_backup;
    Cell[] cell_list_backup;
    size_t cell_count_backup;

    size_t cell_index = CELL_LIST.length;

    Point[] second_prev_row;
     
    while(true) {
	++loop_counter;
	point_list_backup = POINT_LIST;
	point_count_backup = Point.ID_count;
	face_list_backup = FACE_LIST;
	face_count_backup = UFace.face_ID_count;
	cell_list_backup = CELL_LIST;
	cell_count_backup = Cell.cell_ID_count;
	
	try {New_Row = generate_new_row(New_Row);}
	catch(Exception exc){
	    writeln("caught the side node exception");
	    side_node_Flag = 1;}
	npoints_row = New_Row.length;
	writefln("there are %s points in the new paving row", npoints_row);
	assign_node_properties(New_Row);
	adjust_row_ends(New_Row, initial_size, 0);
	point_ID_check();
	

	cell_index = CELL_LIST.length;
	npoints_new_row = New_Row.length;

	//smoothing process:
	biggest_move = smoothing(New_Row, Prev_Row, 1);
	move = biggest_move;
	uint reverse_count = 0;
	for(int i; i<4; ++i){
	    reverse(New_Row);
	    move = smoothing(New_Row, Prev_Row, 1);
	    ++reverse_count;
	}
	if (reverse_count%2 != 0){reverse(New_Row);}

	//smooth all the inside nodes
	//forwards then backwards order to improve symmetry
	inside_nodes = POINT_LIST[npoints_boundary..New_Row[0].ID];
	move = biggest_move;
	reverse(New_Row);
	reverse_count = 1;
	for(int i; i<1; ++i){
	    move = smoothing(inside_nodes, Prev_Row, 0);
	    reverse(New_Row);
	    ++reverse_count;
	}
	if (reverse_count%2 != 0){reverse(New_Row);}
	
	//trialling this function:
	even_out_size(New_Row);
	//has a small effect

	remaining_area = clockwise_check(New_Row); 
	intersections = intersecting_faces(New_Row);
	inside_prev = inside_prev_row_check(New_Row, Prev_Row);

	if(remaining_area > 0 || to!int(intersections.length) > 0 || New_Row.length <= 6 || side_node_Flag == 1 || inside_prev == 0){
	    POINT_LIST = point_list_backup;
	    FACE_LIST = face_list_backup;
	    CELL_LIST = cell_list_backup;
	    Point.ID_count = point_count_backup;
	    UFace.face_ID_count = face_count_backup;
	    Cell.cell_ID_count = cell_count_backup;
	    New_Row = Prev_Row;
	    break;
	}
	//breaks loop at first case of new row intersecting itself
	//or if new row completely inverts
	//or has 6 or less nodes remaining in the paving boundary
	//or if there are no longer 2 consecutive side nodes
	//or if a point in the new row is outside the previous row
	
	//if (loop_counter == 1){break;}
	
	//final smooth
	move = biggest_move;
	reverse_count=0;
	for(int i; i<2; ++i){
	    move = smoothing(New_Row, Prev_Row, 1);
	    reverse(New_Row);
	    ++reverse_count;
	}
	if (reverse_count%2 != 0){reverse(New_Row);}
        
	second_prev_row = Prev_Row;
	Prev_Row = New_Row;
    }

    //trialling this function:
    //even_out_size(New_Row);

    move = biggest_move;
    uint reverse_count=0;
    Point[] prev;
    if (second_prev_row.length != 0) {
	prev = second_prev_row;
    } else { prev = Prev_Row;}
    for(int i; i<2; ++i){
	move = smoothing(New_Row, prev, 1);
	reverse(New_Row);
	++reverse_count;
    }
    if (reverse_count%2 != 0){reverse(New_Row);}

    npoints_row = New_Row.length;
    string reason;
    /* if seam_next = 0 the code proceeds to the closing procedure */
    if(to!int(intersections.length) > 0){
	seam_next = 1; reason = "intersection detected";
	if(inside_prev == 0){
	    reason ~= " because a point in the new row was outside the previous row";}
    }
    else if(remaining_area > 0){
	seam_next = 0; reason = "the new paving row fully inverted";}    
    else if (loop_counter <= 1){
	seam_next = 1; reason="first new row was rejected (no progress)";}
    //not entirely sure if this condition should go here

    if(npoints_new_row <= 6){
	seam_next = 0; reason="there were 6 or less nodes in the new paving row";}
    else if(side_node_Flag == 1){
	seam_next = 0; reason="there are not 2 consecutive side nodes in the paving boundary";}

    writefln("reverted back to %s points in the new paving row because: \n %s", npoints_row, reason);

    assign_node_properties(New_Row);
    adjust_row_ends(New_Row, initial_size, 0);
    point_ID_check();
    Prev_Row = New_Row;
}

void seaming_procedure(ref Point[] New_Row)
{
    Point[][] inner_nodes;
    size_t[] adjust_IDs;
    Point[] revised_New_Row;
    size_t original_nodes_left = New_Row.length;
    double smallest_angle = 2*PI;
    Point smallest_angle_node;
    foreach(node; New_Row){
	if (node.interior_angle < smallest_angle){
	    smallest_angle = node.interior_angle;
	    smallest_angle_node = node;
	}
	if(node.interior_angle < PI/2 && original_nodes_left > 2 ){
	    inner_nodes ~= [seam(node, New_Row, adjust_IDs, original_nodes_left, 0)];
	}
    }
    if (inner_nodes.length == 0){
	writeln("seaming from smallest angle node:");
	inner_nodes ~= [seam(smallest_angle_node, New_Row, adjust_IDs, original_nodes_left, 0)];
    }
   
    if (inner_nodes.length > 0){
	for(size_t i; i<inner_nodes.length; ++i){
	    revised_New_Row ~= New_Row[adjust_IDs[2*i]];
	    revised_New_Row ~= inner_nodes[i];
	    if(i!=inner_nodes.length-1){
		revised_New_Row ~= New_Row[adjust_IDs[2*i+1]..adjust_IDs[2*i+2]];
	    } else {
		size_t row_ID = adjust_IDs[2*i+1];
		while(row_ID != adjust_IDs[0]){
		    revised_New_Row ~= New_Row[row_ID];
		    row_ID = New_Row[row_ID].pointR.rowID;
		}
	    }
	}
	New_Row = revised_New_Row;
        size_t npoints_row = New_Row.length;
	writefln("there are %s points in the new paving row after seaming process", npoints_row);
	//reset row IDs and classifications:
	size_t rowID;
	foreach(node;New_Row){
	    node.store_rowID(rowID);
	    node.node_cat = "none";
	    ++rowID;
	}
     
	assign_node_properties(New_Row);
	adjust_row_ends(New_Row, initial_size, 1);
        
	clockwise_check(New_Row);
    }
}

/*-----------------------------------Section 7----------------------------------*/
/*-----------------------------The Paved Grid Class-----------------------------*/
/*---------------------------------- 150 Lines ---------------------------------*/

class PavedGrid
{
public:
    int dimensions = 2;
    string label;

    size_t nvertices, nfaces, ncells, nboundaries;

    
    this(double[][] boundary_points)
    {
	/*---------------import the boundary geometry---------------------*/
	writeln("importing external boundary...");
	//chose an example geometry file:
	//import_boundary_points("eg_geoms/exterior_surface2.txt");
	foreach(p;boundary_points){
	    make2DPoint(p[0], p[1]);
	}
    	foreach(point;POINT_LIST){
	    point.store_rowID(point.ID);
    	}

	Point[] external_boundary = POINT_LIST;    
	size_t npoints_boundary = external_boundary.length;
	writefln("there are %s points in the external boundary", npoints_boundary);
	if(npoints_boundary%2 != 0){
	    throw new Error(text("imported boundary MUST contain an even number of nodes"));}

	/*----calculate the node information for the exterior boundary----*/
	assign_node_properties(external_boundary);
	adjust_row_ends(external_boundary, initial_size, 1);

	/*----
	  check the imported boundary is: anti-clockwise, closed, non-intersecting
	  ----*/
	if(clockwise_check(external_boundary) > 0) {
	    throw new Error(text("imported external boundary must be anti-clockwise"));
	} else {writeln("boundary is anti-clockwise and closed - good");}

	if(to!int(intersecting_boundary(external_boundary).length) > 0){
	    throw new Error(text("imported boundary intersects itself - not allowed"));
	} else {writeln("boundary does not self-intersect - good");}

	/*--------calculate initial properties for later use------------*/
	initial_size = average_face_length(external_boundary);
	avr_area = initial_size^^2;

	/*---------------generate the first new row---------------------*/
	writeln("generating first new row..."); 
	Point[] New_Row = generate_new_row(external_boundary);
	size_t npoints_row = New_Row.length;
	writefln("there are %s points in the first new paving row", npoints_row);
	assign_node_properties(New_Row);
	Point[] Prev_Row = New_Row;
	adjust_row_ends(New_Row, initial_size, 0);
	smoothing(New_Row, Prev_Row, 1);    
	//this 'ghost' smooth works in the opposite way to the 
	//subsequent paving row smooths
	even_out_size(New_Row);

	/*----
	  continue generating/seaming rows until closing conditions are reached
	  ----*/
	bool seam_next = 1;
	int loop_counter;
	while(true){
	    ++loop_counter;
	    writeln("paving new rows until intersection...");
	    pave_rows_until_intersection(New_Row, Prev_Row, npoints_boundary, seam_next);
	    //if(loop_counter == 1) {break;}
	    if (seam_next==0){break;}
	    writeln("commencing seaming procedure...");
	    seaming_procedure(New_Row);
	    Prev_Row = New_Row;
	    writeln("...completed seam.");
	    //if(loop_counter == 3) {break;}
	}
    	
	/*---------------close the final region-------------------------*/
        
	writeln("commencing closing seam...");
	closing_seam(New_Row);    
	writeln("completing final smoothing...");
	//smooth all the inside nodes:    
	Point[] inside_nodes = POINT_LIST[npoints_boundary..$];
	double biggest_move = smoothing(inside_nodes, Prev_Row, 0);
	double move = biggest_move;
	for(int i; i<2; ++i){
	    move = smoothing(inside_nodes, Prev_Row, 0);
	}
	writeln("DONE.");
        writeln("total cell count: ", CELL_LIST.length);
    	

	/*-------section for size adjustments: not currently included----*/
	/*this works except I think it is somehow generating extra points 
	  which is very difficult to fix, and not sure if excess points 
	  (not utilized by cells) is unacceptable*/

	/*
	  double worst_ratio =  biggest_cell(CELL_LIST, initial_size)/(initial_size^^2);
	  if(worst_ratio > 2){
	  writeln("seam big cells");
	  seam_big_cells(npoints_boundary, initial_size, 0);
	  }

	  writeln("biggest cell ratio: ", worst_ratio);
	  writeln("biggest cell: ", biggest_cell(CELL_LIST, initial_size));
	  writeln("initial_area: ", initial_size^^2);
	*/

	this.nvertices = POINT_LIST.length;
	this.nfaces = FACE_LIST.length;
	this.ncells = CELL_LIST.length;
	this.nboundaries = npoints_boundary;

    } //end: of constructor
} //end: of class definition







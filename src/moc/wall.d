/*
 * file: wall_line.d
 * location: 
 * author: Momar Hughes (adapted from moc_wall.c and moc_wall.h by PA Jacobs)
 * description: set of Wall functions for IMOC_D
 * This is a temporary file for greatly simplified, 2-point linear walls.
 * It will be replaced by wall.d once the complete bezier spline architecture
 * is ready.
 * version:	27 Mar 2015: initial port
 *		10 May 2015: use existing classes in CFCFD3 - Vector3, Line
 *		15 May 2015: simplification to 2-point walls
 */

module wall;

import std.stdio, std.math, std.conv, std.string, std.file;
import kernel;
import geom;
import gpath;
import linesearch;

enum int NWALL = 2; //max number of possible walls
enum int MAX_POINT = 20; //max number of points on each wall
int[NWALL] usable; //flag: 1/0 usable/nonusable Wall
Vector3 WallPoints[NWALL][]; //class from bezier.d
Path[NWALL] WallPath; //class from bezier.d

int WallIsPresent(int iw)
{
    if(iw>=0 && iw<NWALL){
    	return usable[iw];
    } else {
    	return NO;
    } // end if
} // end WallIsPresent()

int WallGetNumberOfPoints(int iw)
{
    // Returns number of points defining the wall AFTI.
    if (WallIsPresent(iw)==YES){
	return to!int(WallPoints[iw].length);
    } else {
	return MOC_ERROR;
    } // end if
} // end WallGetNumberOfPoints()

Vector3 WallGetPoint(int iw,int ip){
    if(WallIsPresent(iw)==YES){
	if(WallPoints[iw][ip] != Vector3()){
	    return WallPoints[iw][ip];
	} else {
	    return Vector3();
	} // end if
    } else {
	return Vector3();
    } // end if
} // end WallGetPoint()

void DeleteWall(int iw){
    if(WallIsPresent(iw)==YES){
	usable[iw]=NO;
	WallPoints[iw] = null;
	WallPath[iw] = null;
    } // end if
} // end DeleteWall()

void WallFromPoints(int iw,string pathstring,Vector3[] points)
{
    if(WallIsPresent(iw)==YES){
	DeleteWall(iw);
    } // end if
    if(points.length>MAX_POINT || points.length<2){
	throw new Error(text("WallFromPoints: invalid number of points: ",points.length));
    } else {
      	if(pathstring=="Line" && points.length==2){
	    WallPath[iw] = new Line(points[0],points[1]);
      	} else if(pathstring=="Arc" && points.length==3){
	    WallPath[iw] = new Arc(points[0],points[1],points[2]);
      	} else if(pathstring=="Arc3" && points.length==3){
	    WallPath[iw] = new Arc3(points[0],points[1],points[2]);
      	} else if(pathstring=="Bezier" && points.length>=2){
	    WallPath[iw] = new Bezier(points);
      	} else if(pathstring=="Polynomial" && points.length>=2){
	    WallPath[iw] = new Polynomial(points);
      	} // end if
    } // end if
    foreach(point;points){
	WallPoints[iw] ~= point;
    } // end foreach*/
    usable[iw] = YES;
} // end WallAddPoint()

void WallFromPolynomial(int iw,double[] coeffs,double x0,double x1)
{
    if(WallIsPresent(iw)==YES){
	DeleteWall(iw);
    } // end if
    Polynomial wall = new Polynomial(coeffs,x0,x1);
    WallPath[iw] = wall;
    WallPoints[iw] ~= wall.P[0];
    WallPoints[iw] ~= wall.P[1];
    usable[iw] = YES;
} // end WallFromPolynomial()

/*
  void WallFromPath(int iw,Path path){
  if(WallIsPresent(iw)==YES){
  DeleteWall(iw);
  } // end if
  auto pathstring = typeof(path).stringof;
  if(pathstring=="Line"){
  Line line = new Line(path);
  WallPoints[iw] = [line.p0,line.p1];
  } else if(pathstring=="Arc"){
  Arc arc = new Arc(path);
  WallPoints[iw] = [arc.a,arc.b,arc.c];
  } else if(pathstring=="Arc3"){
  Arc3 arc3 = new Arc3(path);
  WallPoints[iw] = [arc3.a,arc3.m,arc3.b];
  } else if(pathstring=="Bezier"){
  Bezier bezier = new Bezier(path);
  WallPoints[iw] = bezier.B;
  } else if(pathstring=="Polynomial"){
  Polynomial poly = new Polynomial(path);
  WallPoints[iw] = poly.P;
  } // end if
  WallPath[iw] = path;
  usable[iw] = YES;
  } // end WallAddPoint()
*/
Vector3 WallPos(int iw, double t){
    /* Evaluate the coordinates of the bezier spline for parameter value t
     * Input: 
     iw : specified wall
     t : parametric distance 0 <= t <= 1
     * Output : 
     Returns the coordinates of the wall spline (x,y) or (0,0) on failure AFTI */
    if(WallIsPresent(iw)==YES){
	return WallPath[iw](t);
    } else {
	throw new Error(text("WallPos: invalid wall, iw: ",iw));
    } // end if
} // end WallPos()

double WallSlope(int iw,double t){
    /* Evaluate the slope (dy/dx) of the bezier spline for parameter value t
     * Input: 
     iw : specified wall
     t : parametric distance 0 <= t <= 1
     * Output : 
     Returns the slope of the wall spline or 0 on failure AFTI */
    if(WallIsPresent(iw)==YES){
	Vector3 deriv = WallPath[iw].dpdt(t);
	return deriv.y/deriv.x;
    } else {
	throw new Error(text("WallSlope: invalid wall, iw: ",iw));
    } // end if
} // end WallSlope()

double dist(Vector3 left,Vector3 right){
    return abs(right-left);
} // end dist()

double WallFindT(int iw,Vector3 point,double cos_th,double sin_th,double t_tol=1.0e-6)
{
    if(WallIsPresent(iw)==YES){
	Vector3 xyP = point;
	double t;
	double tL=0.0,tL_new=0.0; Vector3 xyL=WallPos(iw,tL_new); //set left bounds
	double tR=0.0,tR_new=1.0; Vector3 xyR=WallPos(iw,tR_new); //set right bounds
	double alpha = 0.5;
	double dist_tol=t_tol * dist(xyL,xyR);
	double dist_error;
	int count=0;
	do{
	    count++;
	    //writeln("count=",count);
	    tL=tL_new; xyL=WallPos(iw,tL); //update left bounds
	    tR=tR_new; xyR=WallPos(iw,tR); //update right bounds
	    double dx=xyR.x-xyL.x;double dy=xyR.y-xyL.y;
	    double dxPL=xyP.x-xyL.x;double dyPL=xyP.y-xyL.y;
	    double denom=dx*sin_th - dy*cos_th;
	    if(fabs(denom)<=1.0e-10){
		throw new Error(text("WallFindT: Line is parallel to Wall"));
	    } // end if
	    double beta=(dxPL*sin_th-dyPL*cos_th)/denom;
	    Vector3 xyi = Vector3(xyL.x+beta*dx,xyL.y+beta*dy);
	    //
	    t = tL + beta * (tR - tL); 
	    Vector3 xyb = WallPos(iw, t);
	    dist_error = dist(xyb,xyi);
	    //
	    if(beta < 0.0 ) {
		/* Seems to be left of current range; move tL only
		   and reduce the rate at which future ranges are contracted. */
		tL_new = tL - (1.0 + fabs(beta)) * (tR - tL);
		if(tL_new<0.0){tL_new = 0.0;}
		tR_new = tR;
		alpha *= 0.5;
	    } else if(beta > 1.0 ) {
		/* Seems to be right of current range; move tR only
		   and reduce the rate at which future ranges are contracted. */
		tL_new = tL;
		tR_new = tR + (1.0 + beta) * (tR - tL);
		if(tR_new > 1.0){tR_new = 1.0;}
		alpha *= 0.5;
	    } else { // In between the current range; contract the range.
		tL_new = tL + alpha * beta * (tR - tL);
		tR_new = tR - alpha * (1.0 - beta) * (tR - tL);
	    } // end if
	} while(dist_error>dist_tol && count<30 && fabs(tR-tL)>t_tol);
	return t;
    } else {
	throw new Error(text("WallFindT: invalid wall, iw: ",iw));
    } // end if
} // end WallFindT()

/**
 * Finds the closest position along a wall to a point
 * Input: 
 *   iw: index of wall to save
 *   point: coordinates of point
 *   t_tol: tolerance on parametric value when searching for closest point
 * Output:
 *   Returns parametric value of closest wall point
 */
double ClosestWallPos(int iw,Vector3 point,double t_tol=1.0e-4)
{
    if(WallIsPresent(iw)==YES){
	double f(double t){return abs(WallPos(iw,t) - point);}
	double a=0.0,b=1.0;
	minimize!f(a,b,t_tol);
	return 0.5*(a+b);
    } else {
	throw new Error(text("ClosestWallPos: invalid wall, iw: ",iw));
    } // end if
} // end ClosestWallPos()

/**
 * Finds the minimum distance from a point to a wall
 * Input: 
 *   iw: index of wall to save
 *   point: coordinates of point
 *   t_tol: tolerance on parametric value when searching for closest point
 * Output:
 *   Returns minimum distance to wall
 */
double MinimumDistanceToWall(int iw,Vector3 point,double t_tol=1.0e-4)
{
    if(WallIsPresent(iw)==YES){
	double t = ClosestWallPos(iw,point,t_tol);
	return abs(WallPos(iw,t) - point);
    } else {
	throw new Error(text("MinimumDistanceToWall: invalid wall, iw: ",iw));
    } // end if
} // end MinimumDistanceToWall()

/**
 * Checks whether a point can be considered to be on a wall
 * Input: 
 *   iw: index of wall to save
 *   point: coordinates of point to check
 *   tol: maximum distance away from wall that is still considered
 * 		  to be " on " the wall
 * Output:
 *   Returns 1/0 if point is/is not on the wall
 */
int CheckPointOnWall(int iw,Vector3 point,double tol=1.0e-4)
{
    if(WallIsPresent(iw)==YES){ 
	if(MinimumDistanceToWall(iw,point,tol) <= tol){
	    return YES;
	} else {
	    return NO;
	} // end if
    } else {
	throw new Error(text("CheckPointOnWall: invalid wall, iw: ",iw));
    } // end if
}// end CheckPointOnWall()
		
/**
 * Save all important  data from a text file created using SaveNodes
 * Input: 
 *   iw: index of wall to save
 *   FileName: name of file to save to
 */
void SaveWall(int iw,string FileName)
{
    if(WallIsPresent(iw)==YES){
	File fp = File( FileName, "w"); //write only access to FileName
	fp.writefln("%s , no. of points: %d",WallPath[iw].classString(),WallGetNumberOfPoints(iw));
	foreach(Vector3 p;WallPoints[iw]){
	    fp.writefln("%e %e",p.x,p.y);
	} // end foreach
	fp.close();
    } else {
	throw new Error(text("SaveWall: invalid wall, iw: ",iw));
    } // end if
} // end SaveWall()

/**
 * Load all important  data from a text file created using SaveNodes
 * Input: 
 *   iw: index of wall to load to
 *   FileName: name of file to load from
 */
void LoadWall(int iw,string FileName)
{
    if(exists(FileName) == 0){
	throw new Error(text("LoadWall: cannot find file: ",FileName));
    } else if(isFile(FileName) == 0) {
	throw new Error(text("LoadWall: invalid file format, cannot load from file: ",FileName));
    } else {
	if(iw < 0 || iw >= NWALL){
	    throw new Error(text("LoadWall: invalid wall, iw: ",iw));
	} else {
	    DeleteWall(iw); //clear the existing wall
	    File fp = File(FileName, "r");
	    Vector3[] points;
	    string pathstring;
	    while(!fp.eof()){
		auto line = split(fp.readln());
		if(line.length == 0){
		    break;
		} else if(line.length != 2){
		    pathstring = line[0];
		} else {
		    double x = to!double(line[0]);
		    double y = to!double(line[1]);
		    points ~= Vector3(x,y);
		} // end if
	    } // end while
			
	    WallFromPoints(iw,pathstring,points);
	} // end if
    } // end if
} // end LoadWall()	


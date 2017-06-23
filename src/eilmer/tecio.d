/**
 * Interface for tecio library. 
 *
 * libtecio.so can be found at: http://www.tecplot.com/downloads/tecio-library/
 * A small subset of the library is exposed here:
 *   + tecini142
 *   + teczne142
 *   + tecdat142
 *   + tecnode142
 *   + tecend142
 *
 * The D-wrapped functions are named:
 *   + dtecini142
 *   + dteczne142
 *   + dtecdat142
 *   + dtecnode142
 *   + dtecend142
 * 
 * All details of the libtecio can be found in: 
 * http://download.tecplot.com/360/current/360_data_format_guide.pdf
 *
 * Author: Pierpaolo Toniato
 * Date: 2017-05-14
 *
 * History: 2017-06-21 -- some prettifying by RJG
 */


module tecio;

import std.stdio;
import std.string;
import std.conv;
import std.typecons;


string[] zneLst = [
    "ORDERED", "FELINESEG", "FETRIANGLE", "FEQUADRILATERAL", "FETETRAHEDRON",
    "FEBRICK", "FEPOLYGON", "FEPOLYHEDRON"
];

//
// libtecio exposed functions below
//
extern (C) int tecini142(const char* Title, const char* Variables,
			 const char* FName, const char* ScratchDir,
			 int* FileFormat, int* FileType,
			 int* Debug, int* VIsDouble);

extern (C) int tecend142();

extern (C) int teczne142(const char* ZoneTitle, int* ZoneType, int* IMxOrNumPts,
			 int* JMxOrNumElements, int* KMxOrNumFaces, int* ICellMax, int* JCellMax,
			 int* KCellMax, double* SolutionTime, int* StrandId, int* ParentZone,
			 int* IsBlock, int* NumFaceConnections, int* FaceNeighborMode,
			 int* TotalNumFaceNodes, int* NumConnectedBoundaryFaces,
			 int* TotalNumBoundaryConnections, int* PassiveVarList, int* ValueLocation,
			 int* shaveVarFromZone, int* ShareConnectivityFromZone);

//extern (C) int tecdat142(int* N, float* Data, int* IsDouble); // function overloaded to accept both float and doubles
extern (C) int tecdat142(int* N, double* Data, int* IsDouble); //
extern (C) int tecnode142(int* N, int* ConnData);

/***********************************
 * wrapper around tecini142.Initialise output file 
 *
 * Params
 * ------ 
 * title      : dataset title
 * vars       : space-separated list of variables names
 * fname      : filename. Please use plt and szplt for file ending
 * scratchDir : (optional) = Folder were temp files are written. Default current folder
 * fileFmt    : (optional) = 0 plt binary. 1 szplt, (for very big dataset this is recommended)
 * fileType   : (optional) = 0= full solution, 1= grid only, 2= solution only
 * dbg        : (optional) 0 no debug, 1 more output for debugging
 * VIsDouble  : (optional) = Do not use
 *
 * Example dtecini142("IJ Ordered Zones", "X Y P", "ij_ordered.plt")
 * Returns: exit code 0 successful, -1 not successful
 */

int dtecini142(string title, string vars, string fname,
	       string scratchDir=".", int fileFmt=1,
	       int fileType=0, int dbg=0, int VIsDouble=0)
{
    return tecini142(title.toStringz, vars.toStringz,
		     fname.toStringz, scratchDir.toStringz,
		     &fileFmt, &fileType, &dbg, &VIsDouble);
}

/*******************
 * Check whether an  array is null, and return either a pointer or null 
 * Apparently dlang+libtecio do not like passing a pointer of null
 * var = array of ints
 * Returns:
 * pointer to the array or null
 */
int* checknull(int[] var)
{
    int* ptr_var = (var !is null) ? &var[0] : null;
    return ptr_var;
}

/**************************
 * wrapper around teczne142.Initialise zone section (ie each block data)
 *
 * Params: 
 * zoneTitle        : zone name 
 * imaxOrNumPts     : I if structured or Number of nodes if unstructured
 * jmaxOrNumElems   : J if structured or Number of elements if unstructured
 * kmaxOrNumFaces   : (optional) K if structured or Number of faces if special unstructed (see data)
 * solTime          : (optional) Solution time 
 * valLoc           : (optional) = list of 0 or 1 values to indicate if variables are nodal=1 or cell-centered=0
 * 
 * For all the others parameters see tecplot manual
 *
 * Signature do not match teczne142, eliminated some default values to 
 * reduce mess, plus some sane default have been assumed
 * Returns: exit code 0 succesfull, -1 not successfull
 *
 * Example : dteczne142("Zone A",0,2,2,1,360.0,0,[1,1,0])
 */

int dteczne142(string zoneTitle, int zoneType, int imaxOrNumPts,
	       int jmaxOrNumElems, int kmaxOrNumFaces=0,
	       double solTime=0.0, int strandId=0,
	       int[] valLoc=[1,1,0], int parentZone=0,
	       int numFaceConns=0, int faceNbrMode=0,
	       int totalNumBndryConns=1, int[] passiveVarList=null)
{
    int iCellMax, jCellMax, kCellMax;
    int isBlockFmt = 1; // Always use this.
    int totalNumFaceNodes = 1;
    int numConnectedBndryFaces = 1;
    int[] shareVarFromZone = null;
    int shareConnFromZone = 0;

    return  teczne142(zoneTitle.toStringz, &zoneType,
		      &imaxOrNumPts, &jmaxOrNumElems, &kmaxOrNumFaces,
		      &iCellMax, &jCellMax, &kCellMax,
		      &solTime, &strandId, &parentZone, &isBlockFmt,
		      &numFaceConns, &faceNbrMode, &totalNumFaceNodes,
		      &numConnectedBndryFaces, &totalNumBndryConns,
		      checknull(passiveVarList), checknull(valLoc),
		      checknull(shareVarFromZone), &shareConnFromZone);
}

/**************************
 * wrapper around tecdat142.Write zone data
 *
 * Params: 
 * Data = array of values
 * 
 * Function to be called for hwo many variables as necessary
 *
 * Signature do not match tecdat142, other info can be inferred from Data
 * Returns: exit code 0 succesful, -1 not successful
 *
 * Example :=dtecdat142([0.0,0.1,0.2,0.3])
 */

/* Do we ever use floats anymore? Just allow the double version.
int dtecdat142(float[] Data)
{
    int N = to!int(Data.length);
    int IsDouble = 0;
    //Pass by reference of some  array
    float* ptr_Data = &Data[0];
    int code = tecdat142(&N, ptr_Data, &IsDouble);
    return code;
}
*/
int dtecdat142(double[] Data)
{
    int N = to!int(Data.length);
    int isDouble = 1;
    //Pass by reference of some  array
    double* ptrData = &Data[0];
    return tecdat142(&N, ptrData, &isDouble);

}
/**************************
 * wrapper around tecnod142.Write zone connectivity
 * Params: 
 * Data = array of ints
 * Signature do not match tecnode142, other info can be inferred from Data
 * Returns: exit code 0 succesfull, -1 not successfull
 */
int dtecnode142(size_t[] Data)
{
    int[] a;
    foreach (d; Data) a ~= to!int(d); //tecio is really picky
    int* ptrData = &a[0];
    int N = to!int(Data.length);
    return tecnode142(&N, ptrData);
}

/* tecend142 wrapper. Close the file
 *Returns: exit code 0 succesfull, -1 not successfull
*/
int dtecend142()
{
    return tecend142();
}


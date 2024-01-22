/*
  CatalystAdaptor.d
  D port of the CFullExample from Paraview's catalyst 2 catalogue

  TODO List:
   - MPI stuff
 
  Author: Nick Gibbons
*/

import std.stdio;
import std.conv;
import core.stdc.stdlib;
import core.runtime;
import geom.grid.sgrid;
import flowstate;

struct CatalystGrid
{
  uint NumberOfPoints;
  uint NumberOfCells;
  double* Points;
  long* Cells;
}

void InitializeGrid(CatalystGrid* cgrid, StructuredGrid grid)
{

  assert(grid.dimensions==3);

  // create the points -- slowest in the x and fastest in the z directions
  if (cgrid.Points) free(cgrid.Points);

  uint numPoints = to!uint(grid.nvertices);
  cgrid.Points = cast(double*)malloc(3 * double.sizeof * numPoints);

  uint counter = 0;
  foreach(vtx; grid.vertices) {
      cgrid.Points[counter + 0] = vtx.x;
      cgrid.Points[counter + 1] = vtx.y;
      cgrid.Points[counter + 2] = vtx.z;
      counter += 3;
  }

  cgrid.NumberOfPoints = numPoints;

  // create the cell to vertex mapping array
  if (cgrid.Cells) free(cgrid.Cells);

  uint numCells = to!uint(grid.ncells);
  cgrid.Cells = cast(long*)malloc(8 * long.sizeof * numCells);

  counter = 0;
  foreach(k; 0 .. grid.nkv-1){
      foreach(j; 0 .. grid.njv-1){
          foreach(i; 0 .. grid.niv-1){
              size_t[] vtxs = grid.get_vtx_id_list_for_cell(i, j, k);
              foreach(vid; vtxs) {
                  cgrid.Cells[counter] = vid;
                  counter += 1;
              }

          }
      }
  }
  cgrid.NumberOfCells = numCells;
}

void FinalizeGrid(CatalystGrid* grid)
{
  if (grid.Points) free(grid.Points);
  if (grid.Cells) free(grid.Cells);
}

struct Attributes
{
  // A structure for generating and storing point and cell fields.
  // Velocity is stored at the points and pressure is stored
  // for the cells. The current velocity profile is for a
  // shearing flow with U(y,t) = y*t, V = 0 and W = 0.
  // Pressure is constant through the domain.
  double* velx;
  double* vely;
  double* velz;
  double* Pressure;
  CatalystGrid* GridPtr;
}

void InitializeAttributes(Attributes* attributes, CatalystGrid* grid)
{
  attributes.GridPtr = grid;

  uint numCells = attributes.GridPtr.NumberOfCells;

  if (attributes.velx) free(attributes.velx);
  attributes.velx = cast(double*)malloc(double.sizeof * numCells);

  if (attributes.vely) free(attributes.vely);
  attributes.vely = cast(double*)malloc(double.sizeof * numCells);

  if (attributes.velz) free(attributes.velz);
  attributes.velz = cast(double*)malloc(double.sizeof * numCells);

  if (attributes.Pressure) free(attributes.Pressure);
  attributes.Pressure = cast(double*)malloc(double.sizeof * numCells);
}

void FinalizeAttributes(Attributes* attributes)
{
  if (attributes.velx) free(attributes.velx);
  if (attributes.vely) free(attributes.vely);
  if (attributes.velz) free(attributes.velz);
  if (attributes.Pressure) free(attributes.Pressure);
}

void UpdateFields(Attributes* attributes, FlowState[] flowstates)
{
    uint numCells = attributes.GridPtr.NumberOfCells;
    foreach(i; 0 .. numCells) {
        attributes.velx[i]     = flowstates[i].vel.x;
        attributes.vely[i]     = flowstates[i].vel.y;
        attributes.velz[i]     = flowstates[i].vel.z;
        attributes.Pressure[i] = flowstates[i].gas.p;
    }
}

@nogc extern(C) void do_catalyst_initialization();
@nogc extern(C) void do_catalyt_execute(int cycle, double time, CatalystGrid* grid, Attributes* attribs);
@nogc extern(C) void do_catalyt_finalization();

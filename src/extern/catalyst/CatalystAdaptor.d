/*
  CatalystAdaptor.d
  Adapted from a D port of the CFullExample from Paraview's catalyst 2 catalogue

  Author: Nick Gibbons
*/

import std.stdio;
import std.conv;
import core.stdc.stdlib;
import core.runtime;
import geom.grid.sgrid;
import flowstate;

struct CatalystData
{
  uint NumberOfPoints;
  uint NumberOfCells;
  double* Points;
  long* Cells;

  double* velx;
  double* vely;
  double* velz;
  double* Pressure;
}

void InitializeCatalystData(CatalystData* data, StructuredGrid grid)
{

  assert(grid.dimensions==3);

  if (data.Points) free(data.Points);

  uint numPoints = to!uint(grid.nvertices);
  data.Points = cast(double*)malloc(3 * double.sizeof * numPoints);

  uint counter = 0;
  foreach(vtx; grid.vertices) {
      data.Points[counter + 0] = vtx.x;
      data.Points[counter + 1] = vtx.y;
      data.Points[counter + 2] = vtx.z;
      counter += 3;
  }

  data.NumberOfPoints = numPoints;

  // create the cell to vertex mapping array
  if (data.Cells) free(data.Cells);

  uint numCells = to!uint(grid.ncells);
  data.Cells = cast(long*)malloc(8 * long.sizeof * numCells);

  counter = 0;
  foreach(k; 0 .. grid.nkv-1){
      foreach(j; 0 .. grid.njv-1){
          foreach(i; 0 .. grid.niv-1){
              size_t[] vtxs = grid.get_vtx_id_list_for_cell(i, j, k);
              foreach(vid; vtxs) {
                  data.Cells[counter] = vid;
                  counter += 1;
              }

          }
      }
  }
  data.NumberOfCells = numCells;

  if (data.velx) free(data.velx);
  data.velx = cast(double*)malloc(double.sizeof * numCells);

  if (data.vely) free(data.vely);
  data.vely = cast(double*)malloc(double.sizeof * numCells);

  if (data.velz) free(data.velz);
  data.velz = cast(double*)malloc(double.sizeof * numCells);

  if (data.Pressure) free(data.Pressure);
  data.Pressure = cast(double*)malloc(double.sizeof * numCells);
}

void FinalizeCatalystData(CatalystData* data)
{
  if (data.Points) free(data.Points);
  if (data.Cells) free(data.Cells);
  if (data.velx) free(data.velx);
  if (data.vely) free(data.vely);
  if (data.velz) free(data.velz);
  if (data.Pressure) free(data.Pressure);
}


void UpdateCatalystFieldData(CatalystData* data, FlowState[] flowstates)
{
    uint numCells = data.NumberOfCells;
    foreach(i; 0 .. numCells) {
        data.velx[i]     = flowstates[i].vel.x;
        data.vely[i]     = flowstates[i].vel.y;
        data.velz[i]     = flowstates[i].vel.z;
        data.Pressure[i] = flowstates[i].gas.p;
    }
}

@nogc extern(C) void do_catalyst_initialization();
@nogc extern(C) void do_catalyt_execute(int cycle, double time, CatalystData* data);
@nogc extern(C) void do_catalyt_finalization();

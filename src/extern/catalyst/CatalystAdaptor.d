/*
  CatalystAdaptor.d
  Adapted from a D port of the CFullExample from Paraview's catalyst 2 catalogue

  Author: Nick Gibbons
*/

import std.stdio;
import std.conv;
import core.stdc.stdlib;
import core.runtime;
import fluidblock;

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

void InitializeCatalystData(CatalystData* data, FluidBlock[] localFluidBlocks)
{

  if (data.Points) free(data.Points);
  
  uint numPoints = 0;
  foreach(blk; localFluidBlocks) numPoints += to!uint(blk.vertices.length);
  data.Points = cast(double*)malloc(3 * double.sizeof * numPoints);

  uint counter = 0;
  foreach(blk; localFluidBlocks){
      foreach(vtx; blk.vertices) {
          data.Points[counter + 0] = vtx.pos[0].x;
          data.Points[counter + 1] = vtx.pos[0].y;
          data.Points[counter + 2] = vtx.pos[0].z;
          counter += 3;
      }
  }
  data.NumberOfPoints = numPoints;

  // create the cell to vertex mapping array
  if (data.Cells) free(data.Cells);

  uint numCells = 0;
  foreach(blk; localFluidBlocks) numCells += to!uint(blk.cells.length);
  data.Cells = cast(long*)malloc(8 * long.sizeof * numCells);

  uint offset = 0;
  counter = 0;
  foreach(blk; localFluidBlocks) {
      foreach(cell; blk.cells){
          foreach(vtx; cell.vtx){
              data.Cells[counter] = vtx.id + offset;
              counter += 1;
          }
      }
      offset += blk.vertices.length;
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


void UpdateCatalystFieldData(CatalystData* data, FluidBlock[] localFluidBlocks)
{
    uint offset = 0;
    foreach(blk; localFluidBlocks){
        foreach(i; 0 .. blk.ncells) {
            data.velx[offset+i]     = blk.celldata.flowstates[i].vel.x;
            data.vely[offset+i]     = blk.celldata.flowstates[i].vel.y;
            data.velz[offset+i]     = blk.celldata.flowstates[i].vel.z;
            data.Pressure[offset+i] = blk.celldata.flowstates[i].gas.p;
        }
        offset += blk.ncells;
    }
}

@nogc extern(C) void do_catalyst_initialization();
@nogc extern(C) void do_catalyt_execute(int cycle, double time, CatalystData* data);
@nogc extern(C) void do_catalyt_finalization();

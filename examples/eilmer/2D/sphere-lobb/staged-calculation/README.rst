Notes on the Lobb sphere reacting air example
=============================================
:Author: Rowan J. Gollan
:Date: 2016-01-25

This example simulates an experiment performed by Lobb in a
ballistic range. In this experiment, nylon spheres were fired
into air at hypersonic speeds. The shock detachment distance
was measured for a range of flow conditions. The aim of the
experiment was to provide data for validaton of models for
high-temperature air thermochemistry. Since the shock detachment
distance is sensitive to the thermochemical state of the gas
in the shock layer, it was thought that this would be a useful
measurement for chemistry model validation.

This eilmer simulation example is quite advanced: it makes
use of shock-fitting and runs the calculation in stages,
proceeding from coarser grids to finer grids. For a gentler
introduction to blunt body simulation with eilmer, have a look
at the examples called "n90-cylinder" (fixed grid, reacting
nitrogen) and "cylinder-shockfitting" (dynamic shock-fitted grid,
perfect gas).

Configuring the staged calculation
----------------------------------
This calculation proceeds in stages from coarse resolution grids
to finer resolution grids. The idea is to reduce the simulation
time on finer grids since we can start from a converged simulation
on a coarser grid. For those familiar with multi-grid, you might
think of this as a poor-man's one-pass multi-grid approach.
Since we are running the calculation in stages, the user has the
option to configure the details of the stages. Specifically, for
each stage, the user can choose: grid cell count and simulation
length. To do this, edit the tables in 'run-calculation-in-stages.lua'
called "nCells" and "flowTimes". These two tables should be equal in
length. The entries in "nCells" give the number of cells along grid
direction for each of the stage in the calculation. In this grid,
we have used the same number of cells in the 'i' and 'j' directions.
We typically want the number of cells to increase to make use of
a previous coarser solution to start a finer grid calculation.
Note that in the example here, I have defined the numbers of cells
as a mathematical expression. This allowed me to double the value
for ncells along one direction (and so quadruple the total cell
count) for each subsequent stage.

The accompanying table "flowTimes" instructs how long (in simulation
time) to run each calculation. You'll see that these numbers reduce 
because we are hoping to reduce the required time to convergence by
using a previous coarse grid solution as our initial condition.
The "flowTime" is related to the speed of the free stream and a
characteristic length of the body. In this example, one flowtime
is: diameter of sphere / freestream velocity (D/u_inf). The values
in the "flowTimes" table are multipliers on that single flowtime value.
So, for example, the first entry asks for 15 flow times on the
coarsest grid. We typically require much longer on the first grid since
we need time to set up the shock structure from an impulsive initial
condition. On subsequent calculations, the shock structure is in place,
and so a relatively smaller amount of simulation time is required to
get to a converged solution.

How to run this simulation
--------------------------
A Lua script has been provided to coordinate the running of
this calculation in stages. To execute this script, type:
  
  > ./first-time-prep.sh
  > ./run-calculation-in-stages.lua

During execution, the script will place each of stages of the
calculation in a separate subdirectory. The directories are
named 'stage-1', 'stage-2', etc.

Upon successful completion, the script will do a small amount
of post-processing in order to extract the shock detachment
distance. The distances for each stage of the calculation are
collated in the file 'sim-results.txt'.

Special features of the simulation
----------------------------------
I mentioned earlier that this is an advanced example. In this
section, I provide some discussion about some of the special
features that are used in this example. Eilmer provides a
quite sophisticated means of scripting to aid with both
pre-processing and post-processing. This example makes heavy
use of that scripting capability. If some of this seems a bit
complex or confusing, you are encouraged to master some of
the basic examples then return to study this example.

Configuration of prep script via external file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We would like to run a sequence of calculations that are almost
identical but differ only in the number of cells and the length
of simulation time. If we were to do this by hand, we would
copy our input file each time and make modifications to those
two values: number of cells and length of simulation time. We
would keep everything else in the file identical. However, we
want to run this calculation in a "hands-off" scripted manner.
We could teach our simulation script how to write out the entire
prep script and make the required changes as needed. This seems
a bit error prone and difficult to make changes to other aspects
of the simulation. Instead, we "externalise" the parameters we
want to change in an external file. The trick is to get the
preparation script to read this external file as part of setting
up its configuration. The idea then is that this external file
is small and easy for us to create from our simulation script.
Note, we are still getting our simulation script to write out file.
The difference is we only write a small, easily-managed file that
contains just the parameters that change from stage-to-stage,
as opposed to a full-blown preparation script. Let's see this
in a bit more detail.

Our preparation script is called 'lobb.lua' and the external file
is called 'case.lua'. We call this "case" because it contains the
parameters that particular to a given case. Early on in the
preparation script, there is a call to an in-built Lua function::

  dofile('case.lua')

This function takes the contents of the named file ('case.lua')
and executes it Lua code as if it had been typed right in place
where the function is called. Therefore, the external file needs
to contain valid Lua code. In this case, our valid Lua code in
"case.lua" is very simple: it is just a series of declarative 
statements. Here is the contents of a valid 'case.lua' file ::

ncells = 60
no_flow_times = 2.0
oldSoln_jobname = 'lobb'
oldSoln_dir = '../stage-1'
oldSoln_tindx = 100

When I said earlier that we only needed to configure the number
of cells and no_flow_times, that wasn't the whole truth. The extra
information is provided for the staged calculations. For example,
this 'case.lua' file would be used at the second stage of the 
calculation. It instructs the prepration script where to look for
the previous stage (oldSoln). Here it says: "Look for a job named
'lobb' in the 'stage-1' directory, and grab the time index equal
to 100." If that information is *not* present, then the preparation
script assumes we are starting from scratch and will provide a
default inflow boundary to begin the shock-fitting calculation.
For the subsequent calculations, no boundary is required because
we pick up information about the flow and the grid from the old
solution. We discuss this transfer of solution information from
old grid to new in the next section.

To summarise the use of an external configuration file, what
we have done here is "externalise" those aspects of the setup
that change from stage to stage. The rest of the preparation script
is identical for each stage of the preparation. We could write
the external file ('case.lua') by hand. By keeping is simple,
we can also write the file in an automated way. That's what happens
in this example: the master script 'run-calculation-in-stages.lua'
has a function 'prepare_case_file()' that auto-generates a 'case.lua'
file for each particular stage of the calculation.

Transfer of flow solution from coarser simulation to finer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the first stage calculation, we start each subsequent stage
by using the previous coarser grid solution as an initial condition.
As mentioned earlier, this is crude, yet still effective, form of
multi-grid. Eilmer provides two options for setting an initial flow
conditions: 1) a fixed condition set from a `FlowState` object;
or 2) a user-provided function that returns the flow state as a function
of x, y and z. This second option allows for a spatially-varying
initial condition. Since using a coarse grid solution is a
spatially-varying initial condition, we will make use of this second option
in order to initialise our new simulation from an old simulation.

We show here some Lua code that can be used to achieve this solution
transfer as an initial condition. We'll provide some discussion after.
In this code snippet, I have used concrete names for the variables to keep
things simple. Of course, one could parameterise the variable names, and
that has been done in the "lobb.lua" file. ::

   fsol = FlowSolution:new{jobName='lobb', dir='../stage-1',
			   tindx=100, nBlocks=4}
   function initial(x, y, z)
      cell = fsol:find_nearest_cell_centre{x=x, y=y, z=z}
      cell.fmt = "FlowState"
      return fsol:get_cell_data(cell)
   end

What we are interested in here is setting up a function of x, y and z
that returns the flow state each time the function is called. Here that
function is called `initial`. This function depends on a `FlowSolution` object
which is the real workhorse in this procedure.
A `FlowSolution` object is initialised based on a job name, a directory where the
job was run, a particular time index and information about how many blocks in that
old simulation. The `FlowSolution` object will then pick up the flow data at that
particular time index. It makes the data available for inspection and manipulation
by the user in their Lua scripts. The `FlowSolution` object was originally devised
for use in custom post-processing. Here we make use of inspection services it provides.
Let's keep this thought in mind --- that a `FlowSolution` object is initialised
and allows us to inspect its data via some service methods --- as we go on to discuss
the `initial` function.

The `initial` function is what I call a *fillFunction* because it is used to *fill* the
domain with an initial solution. We need to understand a little bit about how that works
and the rules governing a fill function in order to make sense of the `initial` function.
During block initialisation, one of the parameters eilmer needs to know is how to
fill the block with an initial condition. As mentioned earlier, that fill condition could
be a fixed `FlowState` or a it could be a function. If eilmer detects a function, it
then starts looping over all cells in the block (it knows where the cells are based on the
grid that has been provided). Each cell has an (x,y,z) coordinate position. Eilmer then
takes that position and passes it to the user's *fillFunction*. That's why the user's
function must accept three parameters, *x*, *y*, and *z* (and in that order) because
this is what eilmer expects. Eilmer also expects are `FlowState` object or
equivalent table in return.
These are the stipulations that eilmer makes. The user is free to do whatever they want
within their function --- it can be as complex or simple as they like --- so long
as they play by those rules: accept three arguments *x*, *y* and *z*, return a `FlowState`
object or equivalent.
We can see then that our `initial` function fits that pattern. We have defined
`initial` as a function of *x*, *y* and *z*. Then on the first line of the function
we use those arguments and pass them to the inspection method `find_nearest_cell_center{}`.
Based on those (x,y,z) coordinates, `find_nearest_cell_centre{}` then searches the flow
solution to find the nearest cell to that set of (x,y,z) coordinates. Upon locating that cell,
it returns a table that indicates which block and which cell index in that block it found.
The remaining lines are used to massage the data such that we fulfil our obligation of 
returning a `FlowState` object. The return table from `find_nearest_cell_centre{}` contains
a block index *ib* and a cell index *i*. On the following line, we add to that table
(which is a nice feature of Lua) a new field. We add the field *fmt* and give it the value
`FlowState`. This is important for the next step. In that next step, we use another of the
`FlowSolution` inspection methods: `get_cell_data{}`. Having located our cell by block index
and cell index, we pass that information to `get_cell_data{}` so that we can get complete
flow information. That flow information can be returned in various formats. Since we set
the *fmt* type to "FlowState", we get information back as a table that is ready for use
as a `FlowState` object.


Transfer of grid from coarser simulation to finer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On a fixed grid simulation, what we have already discussed about using a *fillFunction*
would be all we need to do in order to transfer an older solution. In this shock-fitting
example, there is an extra complication: the grid has also moved as part of the simulation
process. As such, we also need to build the new grid for a finer resolution case based
on the final grid used on the old solution. We will discuss how to that here.

The key to building a structured grid in eilmer is to provide a parametric surface object.
The parametric surface object is an eilmer provided object that can give (x,y,z) coordinates
as a function of parameters *r* and *s* (in 2D). One of the parametric surface types provided
by eilmer is a `MeshPatch`. The `MeshPatch` accepts a discrete structured grid as input and
builds a mathematical expression that interpolates positions in that grid as a function of
*r* and *s*. Now the interesting bit is that we can use the *old* grid as input to `MeshPatch`.
This ensures that the new grid will have the same edge locations as the old grid. Let's repeat
that one more time because it's a little bit to digest all at once. We need to construct a new
grid that might have more or fewer cells compared to the old grid (so we can't just use the old
grid directly), and we want the new grid to have the same geometry as the old grid. We need to
provide the new grid with a parametric surface. We provide it with a `MeshPatch`. Now a `MeshPatch`
is, in turn, based on a pre-existing structured grid. We give `MeshPatch` the old grid to work with.

Here's how code to do that for a transfer from single-block grid to single-block grid would look.
First we assume we already have our old `FlowSolution` object initialised as *fsol* ::

   oldGrid = fsol:get_sgrid{ib=0}
   psurf = MeshPatch:new{sgrid=oldGrid}
   newGrid = StructuredGrid:new{psurface=psurf, niv=41, njv=41}

Now when construction our new block, we hand it the *newGrid*.

In this example, we have an extra complication to deal with. We have decided to run the simulation
as a multiple block simulation so that we can get some parallel processing benefit of using four
cores of a multi-core workstation. This means we have constructructed four blocks that are stacked
on top of each other in the body tangential direction. It is possible to use our grid transfer technique
shown above directly on each of the four blocks in sequence. However, it gets a bit messy when want to
deal with arbitrary numbers of cell counts. Say the user asks for 90 cells in the tangential direction.
How do we easily split that up over the four blocks? It can be done but it will lead to some mismatch 
in cell sizing at block boundaries. The mismatch will only get worse as the simulation proceeds with
the transfer process for stage after stage. A neater solution is to take the old grid in four pieces,
recombine it as a single grid, do the transfer to the new grid as a single piece, the split the new
grid into four pieces. That way we can make use of our `FluidBlockArray` function which neatly takes
care of splitting 90 cells across four blocks, for example. The code to this join-transfer-split
arrangement is shown here. I don't think it's that much more complicated once you know what's going
on ::

   oldGrid = fsol:get_sgrid{ib=0}
   for i=1,njb-1 do
      print("Joining grid= ", i)
      oldGrid:joinGrid(fsol:get_sgrid{ib=i}, "north")
   end
   psurf = MeshPatch:new{sgrid=oldGrid}
   newGrid = StructuredGrid:new{psurface=psurf, niv=91, njv=91}

The key method in this process is the `joinGrid()` method. What we do on the first line
is initialise the `oldGrid` based on the grid for block 0 of the old solution only. Then
in the for loop, we add the grid from blocks 1, 2, 3 in order to the `oldGrid`. This
adding of grids is achieved with the `joinGrid()` method. Finally, we can create a `MeshPatch`
from the completely joined `oldGrid`.
   
Post-processing to extract shock location
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In these simulations, our metric of interest for comparing to the experimental results is the
distance between the body and the edge of the shock along the stagnation streamline. This is
called the shock detachment or shock standoff distance. To determine this value, on a fixed grid,
we would usually extract the line of data along the stagnation streamline. Then then search along
the line of data from the free stream working towards the body looking for the shock location as
a certain jump in flow values from that of the free stream. In a shock-fitting simulation, as
we have here, the shock location actually comes directly from the calculation and is represented
by the edge of the grid. This actually makes it easier to determine the shock detachment distance.
We need only find the location of the grid vertex at the edge of the grid. To do that, we use a 
custom post-processing script that makes use of the `FlowSolution` object to peek at the data
and grab the vertex value that represents the shock detachment distance. This custom script is 
shown in full here ::

   config.grid_motion = "shock_fitting"
   jobName = "lobb"
   Db = 0.5 * 0.0254 -- diameter (in m) of ball bearing
   -- Pick up flow solution at final time
   fsol = FlowSolution:new{jobName=jobName, dir=".", tindx="last", nBlocks=4}
   vtx = fsol:get_vtx{ib=0, i=0, j=0}
   delta = -vtx.x
   d_D = delta/Db
   f = io.open("shock-detachment.txt", 'w')
   f:write(string.format("%20.12e %20.12e\n", delta, d_D))
   f:close()
   print("shock-detachment= ", delta)
   print("delta/D= ", d_D)

We begin by setting the `grid_motion` configuration option to *shock_fitting*. This ensures that
eilmer is aware that we are working with a moving grid simulation. If we don't do that, eilmer
would assume a fixed grid and would pick up the grid associated with *tindx=0* rather than the
actual grid we need. The main trick in this custom post-processing script are the lines::
   fsol = FlowSolution:new{jobName=jobName, dir=".", tindx="last", nBlocks=4}
   vtx = fsol:get_vtx{ib=0, i=0, j=0}
With these, we pick up the final flow solution and then use the `get_vtx` method to retrieve
the vtx value at the bottom left corner of the grid. The x-ordinate of that location corresponds
to the shock detachment distance. The remainder of the script is used to prettify the output and
convert it to a non-dimensional form (by dividing by the sphere diameter).

We can run this customised script using the `custom-post` option in eilmer::

   e4shared --jobb=lobb --custom-post --script-file=shock-detachment.lua




 
